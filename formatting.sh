#!/usr/bin/env zsh  

# A script to apply clang-format to all source files (*.h, *.cpp) in the folders
# ./src, ./apps, and ./tests.
# Uses gnu_parallel if available. 



gp_command=gnu_parallel
cf_command=clang-format-3.8

echo ""
if [[ $# > 0 ]]; then
    cf_command=$1
    echo -e "\tclang-format is set to ${cf_command}"
fi

if [[ $# > 1 ]]; then
    gp_command=$2
    echo -e "\tgnu_parallel is set to ${gp_command}"
fi

echo ""

if ! type $cf_command > /dev/null; then
    echo -e "\tCould not find '$cf_command'"
    echo -e "\tPlease provide name of clang-format executable as first"
    echo -e "\tpositional argument."
    exit 1
fi

if type $gp_command > /dev/null; then
    echo -e "\tFound gnu_parallel (${gp_command}) ..."
    echo -e "\tFormatting files ... "

    $gp_command $cf_command -style=file -i {}  ::: ./(src|apps|tests)/**/*.h
    $gp_command $cf_command -style=file -i {}  ::: ./(src|apps|tests)/**/*.cpp

else
    echo -e "\tCould not find gnu_parallel (${gp_command}) ..."
    echo -e "\tDoing this the old-fashioned way ..."

    find ./(src|apps|tests) -iname "*.cpp" -or -iname "*.h" -exec $cf_command -style=file -i {} \; 
fi

echo -e "\tDone!"


#gnu_parallel clang-format-3.8 -style=file -i {}  ::: ./(src|apps|tests)/**/*.h
#gnu_parallel clang-format-3.8 -style=file -i {}  ::: ./(src|apps|tests)/**/*.cpp

#gnu_parallel clang-format-3.8 -style=file -i {}  ::: ./apps/**/*.h
#gnu_parallel clang-format-3.8 -style=file -i {}  ::: ./apps/**/*.cpp

