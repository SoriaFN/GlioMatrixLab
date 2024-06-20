//Federico Soria (2024)

dir=getDirectory("Choose the folder with your 2P images");
list=getFileList(dir); 
run("Options...", "iterations=1 count=1 black");

//DIALOG
Dialog.create("Choose your destiny...");
Dialog.addNumber("Acquisition speed in Hz", 1.01);
Dialog.addNumber("Pixel size in um", 0.44);
Dialog.addCheckbox("Eliminate frames at both ends of the stack?", false);
Dialog.addNumber("Number of frames to eliminate", 30);
Dialog.show();
hz=Dialog.getNumber();
pixelWidth=Dialog.getNumber();
cutframe_yes=Dialog.getCheckbox();
cutframe_start=Dialog.getNumber();
frame_duration=1/hz;

setBatchMode(true);

for(i=0;i<list.length;i++){
	filename=dir+list[i];
	if (endsWith(filename, "tif")) {
		open(dir+list[i]);
		name=getTitle();
		setOption("ScaleConversions", true);
		run("16-bit");
		getDimensions(width, height, channels, slices, frames);
		cutframe_end=slices-cutframe_start;
		if (cutframe_yes==true) {
			run("Duplicate...", "duplicate range="+cutframe_start+"-"+cutframe_end);
			print(name+" reduced in "+cutframe_start+" frames at both ends");
		}
		newframes=(cutframe_end-cutframe_start)+1;
		Stack.setXUnit("um");
		Stack.setYUnit("um");
		if (cutframe_yes==true) {
				run("Properties...", "channels=1 slices=1 frames="+newframes+" pixel_width="+pixelWidth+" pixel_height="+pixelWidth+" voxel_depth=1 frame=["+frame_duration+" sec]");
		} else 	run("Properties...", "channels=1 slices=1 frames="+slices+" pixel_width="+pixelWidth+" pixel_height="+pixelWidth+" voxel_depth=1 frame=["+frame_duration+" sec]");
		saveAs("tiff", dir+File.separator+"16bit_"+name);
		name2=getTitle();
		print(name2+" converted and saved in "+dir);
		run("Close All");
	}
}
print("FLAWLESS VICTORY...");