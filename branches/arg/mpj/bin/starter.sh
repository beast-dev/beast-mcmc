#!/bin/sh
#set -x

if [ $# -lt 5 ]; then
  echo "Usage: ./starter.sh -np <np> -dir <dir_of_jar> <name_of_jar>"
  echo "                    [-sport   val (15000)]" 
  echo "                    [-mpjport val (20000)]" 
  echo "                    [-dport   val (10000)]"
  echo "                    [-dev     val (niodev)]" 
  echo "                    [-jvmargs "val"] 
  echo "                    [-appargs "val"] 
  echo "Note: Values in brackets indicate default values for optional params"
  echo "Note: The first five arguments have to be strictly in order"
  echo "Note: 'MPJ_HOME' variable must be set";
  exit 127
fi

np=$2
dir=$4
jar=$5
shift
shift
shift
shift
shift

sport=15000
mpjport=20000 
dport=10000 
dev=niodev  
jvmargs="-Xdebug"
appargs="default_app_args"

while [ "$*" ] ; do
  case $1 in 
   -sport)
     shift
     sport=$1
     shift
     ;;
          
   -mpjport)
     shift
     mpjport=$1
     shift
     ;;
       
   -dport)
     shift
     dport=$1
     shift
     ;;
       
   -dev)
     shift
     dev=$1
     shift
     ;;

   -jvmargs)
     shift
     jvmargs=$1
     shift
     ;;
     
   -appargs)
     shift
     appargs=$1
     shift
     ;;

#   *)
#     appargs="$appargs $1"
#     shift
#     ;;

   esac
done

echo "debugging information -- starts "
echo "np = $np"
echo "dir = $dir"
echo "jar = $jar"
echo "sport = $sport"
echo "mpjport = $mpjport"
echo "dport = $dport"
echo "dev = $dev"
echo "jvmargs = <$jvmargs>"
echo "appargs = <$appargs>"
echo "MPJ_HOME = <$MPJ_HOME>"
echo "debugging information -- ends "

java -jar $MPJ_HOME/lib/starter.jar \
     $np $jar $dir $MPJ_HOME -dev $dev -dport $dport -mpjport $mpjport \
     -sport $sport -args $jvmargs -appargs $appargs
