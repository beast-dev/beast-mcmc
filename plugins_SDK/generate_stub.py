#!/usr/bin/env python
import sys, os, shutil, re
from optparse import OptionParser
PROG_NAME = os.path.basename(os.path.abspath(sys.argv[0]))
VERBOSE = False
def debug(msg):
    if VERBOSE:
        sys.stderr.write('%s: %s\n' % (PROG_NAME, msg))
parser = OptionParser()
parser.add_option('-p', '--package', dest='package', 
                        help='The package that contains the Plugin impl')
parser.add_option('-c', '--class', dest='classname', 
                        help='The name of the of class that implements Plugin')
parser.add_option('-x', '--xmlelement', dest='xml', default='', type='str',
                        help='The name of the xml element associated with the parser')
parser.add_option('-v', '--verbose', dest='verbose', default=False, 
                        action='store_true',
                        help='verbose mode')
                        

(options, args) = parser.parse_args()
if options.verbose:
    VERBOSE = True
if not options.package:
    sys.exit('The package options must be specified')
if not options.classname:
    sys.exit('The class options must be specified')

dir = os.path.dirname(os.path.abspath(sys.argv[0]))
template_par = os.path.join(dir, 'templates')
if not os.path.exists(template_par):
    sys.exit('%s does not exist' % template_par)

sp_package = options.package.split('.')

def copy_with_subst(s, d, rep_dict):
    try:
        source = open(s, 'rU')
    except:
        sys.exit('Could not open %s ' % s)
    try:
        dest = open(d, 'w')
    except:
        sys.exit('Could not open %s ' % d)
    pat_sub_list = [(re.compile(r'@%s@' % k), v) for k, v in rep_dict.iteritems()]
    for line in source:
        for p_el in pat_sub_list:
            line = p_el[1].join(p_el[0].split(line))
        dest.write(line)
    dest.close()
    source.close()
replace_dict = {}
replace_dict['PLUGIN_CLASS'] = options.classname
replace_dict['PLUGIN_PACKAGE'] = sp_package[0]
replace_dict['PLUGIN_SUBPACKAGES'] = '.'.join([''] + sp_package[1:] + [''])
replace_dict['PLUGIN_FULL_PACKAGE'] = options.package
replace_dict['PLUGIN_SRC_PATH'] = '''                <include name="${plugin-package}/**/*.java"/>'''
replace_dict['PLUGIN_CLASS_FILE_PATH'] = '''                <include name="${plugin-package}/**/*.class"/>'''
replace_dict['EXAMPLE_XML_FILE'] = 'test%s.xml' % options.classname
if options.xml:
    xml_el = options.xml
elif options.classname.lower().endswith('plugin'):
    xml_el = options.classname[:len('plugin') - 2]
else:
    xml_el = options.classname

replace_dict['PLUGIN_XML_ELEMENT'] = xml_el

replace_dict['PLUGIN_IMPORTS'] = '''import dr.evomodel.substmodel.NucModelType;
import dr.evomodel.substmodel.FrequencyModel;
import dr.evomodel.substmodel.HKY;
import dr.inference.model.Parameter;
import dr.inference.model.Variable;'''
replace_dict['PLUGIN_XML_SYNTAX_RULES'] = '''new ElementRule(FrequencyModel.FREQUENCIES,
                new XMLSyntaxRule[]{new ElementRule(FrequencyModel.class)}),
        new ElementRule("kappa",
                new XMLSyntaxRule[]{new ElementRule(Variable.class)})'''
replace_dict['PLUGIN_XML_PARSER_STUB'] = '''Variable kappaParam = (Variable) xo.getElementFirstChild("kappa");
    FrequencyModel freqModel = (FrequencyModel) xo.getElementFirstChild(FrequencyModel.FREQUENCIES);

    Logger.getLogger("dr.evomodel").info("Creating THMM substitution model. Initial kappa = " +
            kappaParam.getValue(0));

    return new HKY(kappaParam, freqModel);'''
replace_dict['PLUGIN_PARSER_RETURN_TYPE'] = 'HKY.class'
parser_name = 'DummyModelParser'

debug('\n '.join(['%s = %s' % (k, v) for k, v in replace_dict.iteritems()]))


#sys.exit(str(replace_dict))
dest_dir = options.classname
if os.path.exists(dest_dir):
    sys.exit('%s already exsists' % dest_dir)
dir_list = [dest_dir, 'src'] + sp_package
dest_src_dir = os.path.join(*dir_list)
try:
    os.makedirs(dest_src_dir)
except:
    sys.exit('Could not create %s' % dest_src_dir)
dest_example_dir = os.path.join(dest_dir, 'example')
try:
    os.makedirs(dest_example_dir)
except:
    sys.exit('Could not create %s' % dest_example_dir)

plugin_src = os.path.join(template_par, 'DummyPlugin.java')
dest_src = os.path.join(dest_src_dir, '%s.java' % options.classname)
copy_with_subst(plugin_src, dest_src, replace_dict)

plugin_src = os.path.join(template_par, '%s.java' % parser_name)
dest_src = os.path.join(dest_src_dir, '%sParser.java' % options.classname)
copy_with_subst(plugin_src, dest_src, replace_dict)

plugin_src = os.path.join(template_par, 'build.xml')
dest_src = os.path.join(dest_dir, 'build.xml')
copy_with_subst(plugin_src, dest_src, replace_dict)

plugin_src = os.path.join(template_par, 'beast_sdk.properties.in')
dest_src = os.path.join(dest_dir, 'beast_sdk.properties.in')
copy_with_subst(plugin_src, dest_src, replace_dict)

plugin_src = os.path.join(template_par, 'Dummy.xml')
dest_src = os.path.join(dest_example_dir, 'test%s.xml' % options.classname)
copy_with_subst(plugin_src, dest_src, replace_dict)

