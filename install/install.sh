#Note: must be called with a first argument as the cpipe root
ROOT=$1

function join() {
    local IFS=$1
    shift
    echo "$*"
}

function compile {
# $1 is the directory to compile in
    pushd $1

    if [[ -f  configure.ac ]]; then
        autoconf
    fi
    if [[ -f configure ]]; then
        yes | ./configure
    fi
    if [[ -f Configure ]]; then
        ./Configure -d
    fi
    if [[ -f Makefile ]]; then
        make
    fi

    popd
}

echo `pwd`

for DIRECTORY in $ROOT/tools/*/; do
    compile $DIRECTORY
done

# Add all tool directories and bin folders to PATH
export PATH=`join ':' $ROOT/tools/*/`:`join ':' $ROOT/tools/*/bin/`:$PATH
