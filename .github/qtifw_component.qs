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
            let target_dir = installer.value("TargetDir").replace(/\//g,"\\");
            let add_path = "$LiteralPath = '" + target_dir + "\\bin'; $regPath = 'registry::HKEY_CURRENT_USER\\Environment'; $currDirs = (Get-Item -LiteralPath $regPath).GetValue('Path', '', 'DoNotExpandEnvironmentNames') -split ';' -ne ''; $newValue = ($currDirs + $LiteralPath) -join ';'; Set-ItemProperty -Type ExpandString -LiteralPath $regPath Path $newValue; $dummyName = [guid]::NewGuid().ToString(); [Environment]::SetEnvironmentVariable($dummyName, 'foo', 'User'); [Environment]::SetEnvironmentVariable($dummyName, [NullString]::value, 'User');";
            let remove_path = "$LiteralPath = '" + target_dir + "\\bin'; $regPath = 'registry::HKEY_CURRENT_USER\\Environment'; $currDirs = (Get-Item -LiteralPath $regPath).GetValue('Path', '', 'DoNotExpandEnvironmentNames') -split ';' -ne ''; $newValue = ($currDirs.Split(';') | Where-Object { $_ -ne $LiteralPath }) -join ';'; Set-ItemProperty -Type ExpandString -LiteralPath $regPath Path $newValue; $dummyName = [guid]::NewGuid().ToString(); [Environment]::SetEnvironmentVariable($dummyName, 'foo', 'User'); [Environment]::SetEnvironmentVariable($dummyName, [NullString]::value, 'User');";
            component.addOperation("Execute",
            "powershell.exe",
            "-Command",
            add_path,
            "UNDOEXECUTE",
            "powershell.exe",
            "-Command",
            remove_path);
        }
        else {
            let path_append = ' PATH="$PATH:@TargetDir@/bin"\n';
            if(systemInfo.kernelType == "darwin") { component.addOperation("AppendFile", "@HomeDir@/.bash_profile", "export" + path_append); }
            else { component.addOperation("AppendFile", "@HomeDir@/.bashrc", "export" + path_append); }
            component.addOperation("AppendFile", "@HomeDir@/.cshrc", "setenv" + path_append);
            component.addOperation("AppendFile", "@HomeDir@/.zshrc", "export" + path_append);
            component.addOperation("AppendFile", "@HomeDir@/.profile", "export" + path_append);
         }
    }

    if(systemInfo.kernelType === "winnt") {
        let target_dir = installer.value("TargetDir").replace(/\//g,"\\");
        checkboxForm = component.userInterface("IconCheckBoxForm");
        if(checkboxForm && checkboxForm.iconCheckBox.checked) {
            component.addOperation("CreateShortcut",
            "@TargetDir@/bin/mopac.exe",
            "@DesktopDir@/MOPAC.lnk",
            "iconPath=@TargetDir@/mopac.ico");
        }
        checkboxForm = component.userInterface("FileCheckBoxForm");
        if(checkboxForm && checkboxForm.fileCheckBox.checked) {
            component.addOperation("RegisterFileType",
            "mop",
            target_dir + "\\bin\\mopac.exe \"%1\"",
            "MOPAC input file extension",
            "text/plain",
            target_dir + "\\mopac.ico",
            "ProgId=mopac.mop");
        }
    }
}
