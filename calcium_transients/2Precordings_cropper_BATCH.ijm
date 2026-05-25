//Macro to crop and rescale 2P videos from the Femtonics Setup

dir=getDirectory("Choose the folder with your 2P images");
list=getFileList(dir); 
run("Options...", "iterations=1 count=1 black");

Dialog.create("Choose your destiny...");
Dialog.addNumber("Acquisition speed in Hz", 1.01);
Dialog.show();
hz=Dialog.getNumber();
frame_duration=1/hz;

setBatchMode(true);

for(i=0;i<list.length;i++){
	filename=dir+list[i];
	if (endsWith(filename, "tiff")) {
		open(dir+list[i]);
		name=getTitle();
		getDimensions(width, height, channels, slices, frames);
		getPixelSize(unit, pixelWidth, pixelHeight);
		crop_rescale();
		saveAs("tiff", dir+File.separator+"CROP_"+name);
		name2=getTitle();
		print(name2+" cropped and saved in "+dir);
		run("Close All");
	}
}
print("FLAWLESS VICTORY...");

function crop_rescale() {
	run("Specify...", "width="+height+" height="+height+" x=0 y=0 slice=1");
	run("Crop");
	Stack.setXUnit("um");
	Stack.setYUnit("um");
	run("Properties...", "channels=1 slices=1 frames="+slices+" pixel_width="+pixelWidth+" pixel_height="+pixelWidth+" voxel_depth=1 frame=["+frame_duration+" sec]");
}