function compile {
# $1 is the directory to compile in
    pushd $1

    if [[ -f  configure.ac ]]; then
        autoconf
    fi
    if [[ -f configure* ]]; then
        ./configure
    fi
    if [[ -f Makefile ]]; then
        make
    fi

    popd
}

for DIRECTORY in ./@(data|tools)/*/; do
    compile $DIRECTORY
done
