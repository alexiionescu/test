#!/bin/sh

usage()
{
  echo "Usage: gradle   [ -t | --task TASKNAME] [-c | --clean] [-p | --project]"
  exit 2
}

# : means the option has mandatory param, :: has optional param
# --long are separated by comma ex: --options a:bcd::n --long aaa:,bbb,ccc,ddd::,nnn
PARSED_ARGUMENTS=$(getopt -a -n $0 --options t:c --long task:,clean -- "$@")
if [ $? -ne 0 ]; then usage; fi

CLEAN=false
PROJNAME=$(basename $(pwd))
eval set -- "$PARSED_ARGUMENTS"
while :
do
  case "$1" in
    -t | --task)   TASK="$2"     ; shift 2 ;;
    -p | --project)   PROJNAME="$2"     ; shift 2 ;;
    -c | --clean)  CLEAN=true    ; shift ;;
    # -- means the end of the arguments; drop this, and break out of the while loop
    --) shift; break ;;
    # If invalid options were passed, then getopt should have reported an error,
    # which we checked as VALID_ARGUMENTS when getopt was called...
    *) echo "Unknown option $1"; exit 22; break ;;
  esac
done

EXITCODE=0
INSTALLDIR=${PROJNAME}_INSTALL
echo  ">>> $(date) -> BEGIN Gradle TASK=$TASK CLEAN=$CLEAN <<<"
case $TASK in
    "Make")
        if $CLEAN; then
            make clean
        else
            make
            EXITCODE=$?
            # if [ $EXITCODE -eq 0 ]; then
                # sudo cp -f client/poltysd /usr/local/libexec/bluetooth/
                ## sudo cp -f client/poltys.service /lib/systemd/system/
                # sudo cp -f poltys-utils.sh /usr/local/bin
                # sudo systemctl restart poltys
            # fi
        fi
        ;;
    "Automake")
        if $CLEAN; then
            make distclean
        else
            autoreconf --install
            EXITCODE=$?
            if [ $EXITCODE -eq 0 ]; then
                ./configure
                EXITCODE=$?
            fi
        fi
        ;;
    "AutomakeDebug")
        if $CLEAN; then
            make distclean
        else
            autoreconf --install
            EXITCODE=$?
            if [ $EXITCODE -eq 0 ]; then
                ./configure CFLAGS='-g -Wall' --disable-shared
                EXITCODE=$?
            fi
        fi
        ;;
    "DistributeSrc")
        if $CLEAN; then
            echo "nothing to do"
        else
            make distcheck
            EXITCODE=$?
        fi
        ;;
    "DistributeBin")
        if $CLEAN; then
            rm -rf $INSTALLDIR
            rm -f ${PROJNAME}_binaries.tar.gz
            rm -f ${PROJNAME}_install_files.lst
        else
            make DESTDIR=$INSTALLDIR install
            EXITCODE=$?
            cd $INSTALLDIR
            ## EXTRA FILES
            ## create file list; execute only if not exists; edit and remove unwanted files
            ## if [ ! -f ../${PROJNAME}_files.lst ]; then
                find . -type f -print > ../${PROJNAME}_install_files.lst 
            ## fi
            rm -f ../${PROJNAME}_binaries.tar.gz
            tar zcf ../${PROJNAME}_binaries.tar.gz $(cat ../${PROJNAME}_install_files.lst)
            ## deploy on target 
            ## sudo tar -xzf ${PROJNAME}_binaries.tar.gz -C /
            ## sudo systemctl enable myservice
            ## sudo systemctl start myservice
        fi
        ;;
esac

echo ">>> $(date) -> END Gradle TASK=$TASK CLEAN=$CLEAN EXITCODE=$EXITCODE<<<"

exit $EXITCODE