function Component()
{
    installer.installationFinished.connect(this, Component.prototype.installationFinishedPageIsShown);
}

Component.prototype.installationFinishedPageIsShown = function()
{
    try {
        if(installer.isInstaller() && (installer.status == QInstaller.Success)) {
            if(systemInfo.kernelType === "winnt") {
                installer.addWizardPageItem(component, "FileCheckBoxForm", QInstaller.InstallationFinished, 2);
                installer.addWizardPageItem(component, "IconCheckBoxForm", QInstaller.InstallationFinished, 3);
            }
            else {
                installer.addWizardPageItem(component, "PathCheckBoxForm", QInstaller.InstallationFinished, 1);
            }
        }
    }
    catch(e) {
        console.log(e);
    }
}

Component.prototype.createOperations = function()
{
    component.createOperations();

    var checkboxForm = component.userInterface("PathCheckBoxForm");
    if(checkboxForm && checkboxForm.pathCheckBox.checked) {
        if(systemInfo.kernelType === "winnt") {
            component.addOperation("Execute",
            "@TargetDir@/add-to-path.bat",
            "@TargetDir@\bin",
            "UNDOEXECUTE",
            "@TargetDir@/remove-from-path.bat",
            "@TargetDir@\bin");
        }
        else {
            component.addOperation("Execute",
            "bash",
            "@TargetDir@/add-to-path.sh",
            "@TargetDir@/bin",
            "UNDOEXECUTE",
            "bash",
            "@TargetDir@/remove-from-path.sh",
            "@TargetDir@/bin");
            component.addOperation("Execute",
            "bash",
            "@TargetDir@/add-to-path.sh",
            "@TargetDir@/lib",
            "UNDOEXECUTE",
            "bash",
            "@TargetDir@/remove-from-path.sh",
            "@TargetDir@/lib");
         }
    }

    if(systemInfo.kernelType === "winnt") {
        checkboxForm = component.userInterface("IconCheckBoxForm");
        if(checkboxForm && checkboxForm.iconCheckBox.checked) {
            component.addOperation("CreateShortcut",
            "@TargetDir@/bin/mopac-win.exe",
            "@DesktopDir@/MOPAC.lnk",
            "iconPath=@TargetDir@/mopac.ico");
        }
        checkboxForm = component.userInterface("FileCheckBoxForm");
        if(checkboxForm && checkboxForm.fileCheckBox.checked) {
            component.addOperation("RegisterFileType",
            ".mop",
            "@TargetDir@/bin/mopac-win.exe '%1'",
            "MOPAC input file extension",
            "text/plain",
            "@TargetDir@/mopac.ico");
        }
    }
}
