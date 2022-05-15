function Component()
{
    component.loaded.connect(this, addCheckBoxes);
}

addCheckBoxes = function()
{
    if(installer.isInstaller()) {
        installer.addWizardPageItem(component, "PathCheckBoxForm", QInstaller.TargetDirectory, 1);
        if(systemInfo.kernelType === "winnt") {
            installer.addWizardPageItem(component, "FileCheckBoxForm", QInstaller.TargetDirectory, 2);
            installer.addWizardPageItem(component, "IconCheckBoxForm", QInstaller.TargetDirectory, 3);
        }
    }
}

Component.prototype.createOperations = function()
{
    component.createOperations();

    var checkboxForm = component.userInterface("PathCheckBoxForm");
    if(checkboxForm && checkboxForm.pathCheckBox.checked) {
        if(systemInfo.kernelType === "winnt") {
            component.addOperation("Execute",
            "pwsh.exe",
            "@TargetDir@/add-path.ps1",
            "UNDOEXECUTE",
            "pwsh.exe",
            "@TargetDir@/remove-path.ps1");
        }
        else {
            let path_append = ' PATH="$PATH:@TargetDir@/bin"\n';
            component.addOperation("AppendFile", "@HomeDir@/.bashrc", "export" + path_append);
            component.addOperation("AppendFile", "@HomeDir@/.cshrc", "setenv" + path_append);
            component.addOperation("AppendFile", "@HomeDir@/.zshrc", "export" + path_append);
            component.addOperation("AppendFile", "@HomeDir@/.profile", "export" + path_append);
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
