/*
 * Macro template to process multiple images in a folder
 */

#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix

// See also Process_Folder.py for a version of this code
// in the Python scripting language.

processFolder(input);

// function to scan folders/subfolders/files to find files with correct suffix
function processFolder(input) {
	list = getFileList(input);
	list = Array.sort(list);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], suffix))
			processFile(input, output, list[i]);
	}
}

// hop1 analysis
function processFile(input, output, file) {
	path = input + File.separator;
	full_path = path + file; 
	open(full_path); // image is now open

	// set up variables
	title = getTitle();
	

	// select RGB file and split channels for separate analysis
	run("Split Channels");

	// set measurements
	run("Set Measurements...", "area mean standard min integrated display redirect=None decimal=3");
	
	// set each file to 16-bit for consistency (after split channels)

	run("ROI Manager...");
	
	// wait for user input -> select circle around dapi
	// if you don't select anything and press OK there will be an error
	// ** fix this w a loop? 
	
	setTool("freehand"); // select freehand tool
	waitForUser("Draw outline around DAPI, then hit OK"); 

	roiManager("add"); // add selection to ROI Manager
	
	color = newArray("(blue)","(red)","(green)"); // array of strings for individual RGB images
	
	// measure
	for (i=0; i<3; i++) {
		
		selectWindow(title + " " + color[i]);
		//print(color[i]);
		run("16-bit");
		
		roiManager("Select", 0);
		
		run("Measure");

	}
	
	roiManager("Deselect");
	roiManager("Delete");


	dir = getDirectory("home");
	close("*");
}


selectWindow("Results");
saveAs("results",  output + File.separator + "Zip1TC Intensity Data.csv");
close("Zip1TC Intensity Data.csv");
close("Results");
close("ROI Manager");