function Controller()
{
}

Controller.prototype.IntroductionPageCallback = function()
{
    var widget = gui.currentPageWidget(); // get the current wizard page
    if (widget != null && installer.isInstaller()) {
        widget.MessageLabel.setText("Welcome to the MOPAC Setup Wizard.\n\n\
This is the first official release of the open-source version of MOPAC [https://github.com/openmopac/mopac], \
which continues from the commercial development of MOPAC that has ended with MOPAC2016. \
The original MOPAC developer, James J. P. Stewart, will continue to be involved with MOPAC at a gradually diminishing level of activity.\n\n\
For long-time users, there are some small changes to be aware of. \
The main executable has been renamed from \"MOPAC2016\" to \"mopac\", and it now depends \
on a core shared library, \"libmopac\", as well as Intel's OpenMP shared library, \
both of which are contained in this installer package."); // set the welcome text

        installer.setDefaultPageVisible(QInstaller.StartMenuSelection, false);
    }
}
