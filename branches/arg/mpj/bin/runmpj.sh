#!/bin/sh
#set -x

if [ $# -ne 2 ]; then
  echo "Usage: runmpj.sh <conf_file> <jar_file>[OR]<class_file>";
  echo "       For e.g. runmpj.sh ../conf/mpj2.conf ../lib/test.jar";
  echo "       For e.g. runmpj.sh mpj.conf hello.World";
  exit 127
fi

conf=$1
lines=`cat $conf  | egrep -v "#" | egrep "@"`
dir=`pwd`
name=$2
count=0

backslash2slash() {
    echo $1 | sed 's/\\/\//g'
}

CLASSPATH_SEPARATOR=":"
case "`uname`" in
  CYGWIN*) 
    MPJ_HOME=`backslash2slash $MPJ_HOME`
    CLASSPATH=`backslash2slash $CLASSPATH`

    CLASSPATH_SEPARATOR=";"
    ;;
esac


for i in `echo $lines`; do 

  host=`echo $i | cut -d "@" -f 1`
  rank=`echo $i | cut -d "@" -f 3`    

  case "$name" in
    *.jar )

      ssh $host "cd $dir; export MPJ_HOME=$MPJ_HOME ; \
          java -cp \"$MPJ_HOME/lib/mpj.jar${CLASSPATH_SEPARATOR}$CLASSPATH\" \
               -jar $name $count $conf niodev;" &

      ;;

    * )

echo $MPJ_HOME
      ssh $host "cd $dir; export MPJ_HOME=$MPJ_HOME ; \
          java -cp \"$MPJ_HOME/lib/mpj.jar${CLASSPATH_SEPARATOR}$CLASSPATH\" \
               $name $count $conf niodev;" & 
      ;;
  esac

  count=`expr $count + 1`

done
