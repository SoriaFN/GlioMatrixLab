//Initialization
print("\\Clear");
run("Options...", "iterations=1 count=1 black");
roiManager("reset");
run("Clear Results");
run("Set Measurements...", "area center limit display redirect=None decimal=2");
name=getTitle();
dir=getDirectory("image");
waitForUser("Arrange the windows to your liking");

print("Analyzing "+name);
print("Quantifying Total Image Area...");
run("Select All");
run("Measure");
run("Select None");
run("Split Channels");

//Iba1 channel
print("Creating Iba1 mask...");
selectWindow("C1-"+name);
run("Z Project...", "projection=[Max Intensity]");
run("Gaussian Blur...", "sigma=1");
setAutoThreshold("Default dark");
run("Threshold...");
waitForUser("Set Threshold", "Set Threshold level using the UPPER sliding bar \nThen click OK. \n \nDo not press Apply!");
run("Convert to Mask");
run("Median...", "radius=2");
saveAs("tiff", dir+File.separator+"Iba1_BIN_"+name);
print("Iba1 mask saved in "+dir);
print("Quantifying Iba1 Area...");
rename("Iba1_"+name);
setAutoThreshold("Default dark");
run("Measure");
print("Creating selection from Iba1 mask (wait)...");
run("Create Selection");
roiManager("Add");
run("Select None");
selectWindow("Results");
saveAs("Results", dir+File.separator+"Iba1_"+name+".xls");
run("Clear Results");
selectWindow("Iba1_"+name);

//NM channel
print("Creating NM mask...");
selectWindow("C4-"+name);
run("Z Project...", "projection=[Sum Slices]");
selectWindow("SUM_C4-"+name);
run("Despeckle");
saveAs("tiff", dir+File.separator+"NM_"+name);
run("Duplicate...", " ");
setAutoThreshold("RenyiEntropy");
run("Threshold...");
waitForUser("Set Threshold", "Set Threshold level using the LOWER sliding bar \nThen click OK. \n \nDo not press Apply!");
run("Convert to Mask");
run("Median...", "radius=2");
saveAs("tiff", dir+File.separator+"NM_FULL_"+name);
print("NM mask saved in "+dir);
print("Quantifying Total NM (NM_FULL)...");
rename("NM_FULL_"+name);
setAutoThreshold("Default dark");
run("Analyze Particles...", "  show=Nothing display exclude summarize");

//NM within Iba1
print("Quantifying NM within Iba1 mask (NM_IN)...");
selectWindow("NM_FULL_"+name);
run("Duplicate...", " ");
rename("NM_IN_"+name);
roiManager("Select", 0);
run("Enlarge...", "enlarge=1");
run("Clear Outside");
run("Select None");
run("Analyze Particles...", "  show=Nothing display exclude summarize");
saveAs("tiff", dir+File.separator+"NM_IN_"+name);

//NM outside Iba1
print("Quantifying NM outside Iba1 mask (NM_OUT)...");
selectWindow("NM_FULL_"+name);
run("Duplicate...", " ");
rename("NM_OUT_"+name);
roiManager("Select", 0);
run("Enlarge...", "enlarge=0.5");
run("Clear");
run("Select None");
run("Analyze Particles...", "  show=Nothing display exclude summarize");
saveAs("tiff", dir+File.separator+"NM_OUT_"+name);

//MergeImage
selectWindow("NM_"+name);
run("Add Image...", "image=[Iba1_"+name+"] x=0 y=0 opacity=50");
saveAs("tiff", dir+File.separator+"MERGE_"+name);

//Save files
print("Saving Files...");
selectWindow("Results");
saveAs("Results", dir+File.separator+"NM_ALL_"+name+".xls");
selectWindow("Summary");
saveAs("Results", dir+File.separator+"NM_SUMMARY_"+name+".xls");
run("Tile");
print("DONE!!!");
showMessage("FLAWLESS VICTORY");
