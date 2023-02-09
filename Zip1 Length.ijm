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

// dapi analysis
function processFile(input, output, file) {
	path = input + File.separator;
	full_path = path + file; 
	open(full_path); // image is now open

	// set up variables
	title = getTitle();
	title1 = File.nameWithoutExtension;

	// select RGB file and split channels for separate analysis
	run("Split Channels");

	// zip1 analysis
	close(title + " (green)"); // close hop1 channel
	close(title + " (blue)"); // close dapi channel
	selectWindow(title + " (red)"); // open zip1 channel
	run("Auto Threshold", "method=MaxEntropy dark"); //Max Entropy looks most representative for zip1 channel
													 //it doesn't capture everything, but it's the only one that
													 //separates the zip1 tracks enough for there to be individual tracks
	run("Convert to Mask");
	run("Skeletonize");
	run("Analyze Skeleton (2D/3D)", "prune=[shortest branch] show");
	
	dir = getDirectory("home");

	selectWindow(title + " (red)");
	saveAs("Tiff", output + File.separator + "skeleton" + File.separator + title1 + " zip1skeleton");
	close(title1 + " zip1skeleton.tif");
	close("Tagged skeleton");

	close("Results");
	selectWindow("Branch information");
	saveAs("results",  output + File.separator + "branch info" + File.separator + title1 + " Zip1Length.csv");
	close(title1 + " Zip1Length.csv");
	
}



