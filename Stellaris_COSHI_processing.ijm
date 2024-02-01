//GUI
Dialog.create("Reverse option");
Dialog.addCheckbox("Reverse stack?", false);
Dialog.addNumber("Value for Mean Filter", 1);
Dialog.show();
reverse=Dialog.getCheckbox();
value_for_meanfilter=Dialog.getNumber();

//FILE
name=getTitle();
run("Duplicate...", "title=TEMP duplicate");
save_file=getBoolean("Save image?");
if (save_file==true) {
	dir=getDirectory("Choose a folder to save the files");
}

//PROCESSING
if (reverse==true) {
	selectImage("TEMP");
	run("Reverse");
	run("Bleach Correction", "correction=[Histogram Matching]");
	print("Image has been reversed for bleach correction.");
} else {
	selectImage("TEMP");
	run("Bleach Correction", "correction=[Histogram Matching]");
}
waitForUser("Is bleach correction OK?");
selectImage("DUP_TEMP");
if (reverse==true) {
	run("Reverse");
}
run("Subtract Background...", "rolling=500 stack");
run("StackReg ", "transformation=[Rigid Body]");
run("Mean...", "radius="+value_for_meanfilter+" stack");

if (save_file==true) {
	selectImage("DUP_TEMP");
	saveAs("Tiff", dir+File.separator+"REG_SB_BC_"+name);
	print("Processed image saved in "+dir);
}
close_images=getBoolean("Close all images?");
if (close_images==true) {
	run("Close All");
}
print("DONE");
