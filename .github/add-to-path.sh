#!/bin/bash
# add first command line argument to PATH using initialization files

path_append="PATH=\"\$PATH:$1\""

# .bashrc (bash shell)
echo "export $path_append" >> ~/.bashrc

# .cshrc (csh shell)
echo "setenv $path_append" >> ~/.cshrc

# .zshrc (zsh shell)
echo "export $path_append" >> ~/.zshrc

# .profile (other shells)
echo "export $path_append" >> ~/.profile
