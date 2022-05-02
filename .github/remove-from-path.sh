#!/bin/bash
# remove first command line argument from initialization files that add it to the PATH

path_append="PATH=\"\$PATH:$1\""

# .bashrc (bash shell)
echo "$(grep -v "export $path_append" ~/.bashrc)" > ~/.bashrc

# .cshrc (csh shell)
echo "$(grep -v "setenv $path_append" ~/.cshrc)" > ~/.cshrc

# .zshrc (zsh shell)
echo "$(grep -v "export $path_append" ~/.zshrc)" > ~/.zshrc

# .profile (other shells)
echo "$(grep -v "export $path_append" ~/.profile)" > ~/.profile
