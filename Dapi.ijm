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

	// select RGB file and split channels for separate analysis
	run("Split Channels");

	// dapi analysis
	close(title + " (green)"); // close hop1 channel
	close(title + " (red)"); // close zip1 channel
	selectWindow(title + " (blue)"); // open dapi channel
	setAutoThreshold("RenyiEntropy dark");
	run("Convert to Mask");
	run("Analyze Particles...", "display summarize add");
	dir = getDirectory("home");
	close(title + " (blue)");
}

close("Results");
selectWindow("Summary");
saveAs("results",  output + File.separator + "DapiArea.csv");
close("DapiArea.csv");


