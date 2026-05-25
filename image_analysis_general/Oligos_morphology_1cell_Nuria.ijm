/* MACRO NURIA OLIGOS
 * ------------------
 * 
 * Gliomatrix Lab 2024 - federico.soria@achucarro.org
 */

var ch=3;          //This is the oligo's channel
var savefile=true; //Set this to 'false' if you don't want to save files.

//CLEAR PREVIOUS RESULTS, ROI AND LOG INFO
run("Collect Garbage");
roiManager("reset");
roiManager("Show None");
run("Clear Results");
run("Options...", "iterations=1 count=1 black");
setOption("BlackBackground", true);
setForegroundColor(255, 255, 255);
setBackgroundColor(0, 0, 0);

//INITIALIZATION
if (nImages==0) {
	exit("No images open. Please open an image");
}
name=getTitle();
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);
if (channels==1) {
	exit("Image is not multichannel. Please open a multichannel image");
}
if (slices==1) {
	exit("Image is not z-stack. Please open a multislice image");
}
if (savefile==true) {
	dir=getDirectory("Choose your destiny...");
}
run("Set Measurements...", "area fit display redirect=None decimal=2");

//PREPROCESSING
run("Split Channels");
selectWindow("C3-"+name);
run("Z Project...", "projection=[Max Intensity]");
setTool("polygon");
waitForUser("Draw a polygon tightly around the cell. Then, press OK");
run("Clear Outside");
run("Select None");
if (savefile==true) {
	selectWindow("MAX_C3-"+name);
	saveAs("tiff", dir+File.separator+"CROP_C3-"+name);
}
name_max=getTitle();

//SEGMENTATION
selectWindow(name_max);
setAutoThreshold("Default dark");
run("Threshold...");
waitForUser("Set Threshold", "Set Threshold level using the upper sliding bar \nThen click OK. \n \nDo not press Apply!");
run("Convert to Mask");
run("Median...", "radius=1");
run("Analyze Particles...", "size=100-infinity add");
fix=getBoolean("Do you want to manually fix the image?", "YES, Fix it", "NO, continue");
if (fix==true) {
	do{
		roiManager("reset");
		roiManager("Show None");
		run("Invert LUT");
		setTool("Paintbrush Tool");
		waitForUser("Fix the processes with 3px brush.\nUse original image as guide.\n \nALT+CLICK to paint.\nCLICK to clear.");
		run("Invert LUT");
		run("Convert to Mask");
		run("Analyze Particles...", "size=100-infinity add");
		
		cont=getBoolean("Do you want to continue fixing?", "YES, needs more tuning", "NO, it's OK now. Continue");
	}while(cont==true);
} print("Binary image fixed manually.");
selectWindow(name_max);
if (savefile==true) {
	selectWindow(name_max);
	saveAs("tiff", dir+File.separator+"BIN_C3-"+name);
}
name_bin=getTitle();

//QUANTIFICATION OF ASPECT RATIO
roiManager("reset");
run("Analyze Particles...", "size=100-infinity add");
n_cells=roiManager("count");
print(n_cells+" ROIs detected");
if (n_cells>1) {
	waitForUser("WARNING:\n \nBinary image contains more than one ROI.\nOnly first ROI will be taken into account");
}
roiManager("Select", 0);
run("Clear Outside");
run("Convex Hull");
run("Measure");
area=getResult("Area", 0);
major_axis=getResult("Major", 0);
minor_axis=getResult("Minor", 0);
aspect_ratio=minor_axis/major_axis;

//QUANTIFICATION OF SKELETON
selectWindow(name_bin);
run("Duplicate...", "title=SKELETON");
run("Select None");
run("Skeletonize");
run("Analyze Skeleton (2D/3D)", "prune=none calculate");
Table.sort("Longest Shortest Path");
n=Table.size-1;
lsp=getResult("Longest Shortest Path", n);
if (savefile==true) {
	selectWindow(name_bin);
	saveAs("tiff", dir+File.separator+"SKL_C3-"+name);
}

//CUSTOM TABLE
myTable(name,unit,area,major_axis,minor_axis,aspect_ratio,lsp);

function myTable(a,b,c,d,e,f,g){
	title1="Quantification";
	title2="["+title1+"]";
	if (isOpen(title1)){
   		print(title2, a+"\t"+b+"\t"+c+"\t"+d+"\t"+e+"\t"+f+"\t"+g);
	}
	else{
   		run("Table...", "name="+title2+" width=1000 height=300");
   		print(title2, "\\Headings:File\tUnit\tArea\tMajor Axis\tMinor Axis\tAspect Ratio\tLongest Shortest Path");
   		print(title2, a+"\t"+b+"\t"+c+"\t"+d+"\t"+e+"\t"+f+"\t"+g);
	}
}
waitForUser("Copy Results to Excel");

//EXIT
close_images = getBoolean("Close all images?");
if (close_images==true) {
	run("Close All");	
}
if (savefile==true) {
	selectWindow("Log");
	saveAs("Text", dir+File.separator+"LOG_"+name+".txt");
}

print("Log and images saved in "+dir);
print("Quantification Table needs to be manually saved");

