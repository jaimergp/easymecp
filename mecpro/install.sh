#!/bin/bash

#Copyright (c) 2017 Brigham Young University

#See the file license.txt for copying permission.

echo "Creating template..."
if test ! -e ~/.templates; then
	mkdir ~/.templates
fi
sed 's|~HOME~|'"$HOME"'|g' mecpraw.template > ~/.templates/mecp.template

echo "Installing scripts..."
if test ! -e ~/scripts; then
	mkdir ~/scripts
fi

if test ! -e ~/bin; then
	mkdir ~/bin
fi

if test ! -e ~/pylib; then
	mkdir ~/pylib
fi


cp -rf mecpnext.sh ~/scripts/mecpnext
cp -rf mecpstatus.sh ~/scripts/mecpstatus
cp -rf mecpcheck.sh ~/scripts/mecpcheck
cp -rf mecpstart.sh ~/scripts/mecpstart
cp -rf mecpstatus_helper.sh ~/scripts/mecpstatus_helper
cp -rf mecpdata.sh ~/scripts/mecpdata
cp -rf mecpinput.py ~/scripts/mecpinput.py
cp -rf mecp.py ~/bin/mecpopt
cp -rf getwalltime.py ~/scripts/hours2time
cp -rf atom.py gaussian.py rotate.py mecpinput.py ~/pylib/
chmod +x ~/scripts/* ~/bin/*

if ! test -e ~/.bash_profile; then
	echo "#!/bin/bash" > ~/.bash_profile
fi

if ! test -e ~/.modules; then
	echo "#%Module" > ~/.modules
fi

if ! grep -q 'export PATH="$PATH:~/scripts"' ~/.bash_profile; then
	echo "Setting up path..."
	echo >> ~/.bash_profile
	echo 'export PATH="$PATH:~/scripts"' >> ~/.bash_profile
	#This is why you should run this with the source command
	export PATH="$PATH:~/scripts"
fi

if ! grep -q 'export PATH="$PATH:~/bin"' ~/.bash_profile; then
	echo "Setting up path..."
	echo >> ~/.bash_profile
	echo 'export PATH="$PATH:~/bin"' >> ~/.bash_profile
	#This is why you should run this with the source command
	export PATH="$PATH:~/bin"
fi

if ! grep -q 'export PYTHONPATH="$PYTHONPATH:$HOME/pylib"' ~/.bash_profile; then
	echo "Setting up python path..."
	echo >> ~/.bash_profile
	echo 'export PYTHONPATH="$PYTHONPATH:$HOME/pylib"' >> ~/.bash_profile
	#This is why you should run this with the source command
	export PYTHONPATH="$PYTHONPATH:$HOME/pylib"
fi

if ! grep -q "module add python/2.7.5" ~/.modules; then
	echo "Using Python 2.7.5"
	echo >> ~/.modules
	echo 'module add python/2.7.5' >> ~/.modules
	module add python/2.7.5
fi

echo "Installation Complete."

echo "If you didn't run this with the 'source' command, run the following command:"
echo " > source ~/.bash_profile"
echo "Or log out and back in. (whichever you prefer)"
