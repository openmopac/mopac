function Controller()
{
}

Controller.prototype.IntroductionPageCallback = function()
{
    var widget = gui.currentPageWidget(); // get the current wizard page
    if (widget != null && installer.isInstaller()) {
        widget.MessageLabel.setText("Welcome to the MOPAC Setup Wizard.\n\n\
This is the first official release of the open-source version of MOPAC [https://github.com/openmopac/mopac]. \
It is a continuation of the commercial development of MOPAC, which has ended with MOPAC2016. \
The original MOPAC developer, James J. P. Stewart, will continue to be involved with the project at a gradually diminishing level of activity.\n\n\
For long-time users, there are some small but important changes to the distribution of MOPAC. \
While it will continue to be distributed as a stand-alone package, now with a graphical installer, it will also be increasingly available over various package managers. \
Also, the main executable has be renamed to \"mopac\", and it is now dependent on a core shared library, \"libmopac\", as well as Intel's OpenMP shared library, \
both of which are contained in this installer package."); // set the welcome text
    }
}
