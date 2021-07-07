#!/bin/bash

#cp old_newmannziff.cpp newmannziff.cpp
sed -ie'' 's/5/=/g' newmannziff.cpp
sed -ie'' 's/\%/}/g' newmannziff.cpp
sed -ie'' 's/,/</g' newmannziff.cpp
sed -ie'' 's/\./>/g' newmannziff.cpp
sed -ie'' 's/\~/(/g' newmannziff.cpp
sed -ie'' 's/\!/)/g' newmannziff.cpp
sed -ie'' 's/@/[/g' newmannziff.cpp
sed -ie'' 's/\#/]/g' newmannziff.cpp
sed -ie'' 's/\$/{/g' newmannziff.cpp
