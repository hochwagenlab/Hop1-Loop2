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
	title0 = File.nameWithoutExtension + "_Hop1NoThreshold";
	title1 = File.nameWithoutExtension + "_Hop1Threshold";
	title2 = File.nameWithoutExtension + "_Hop1Foci.csv";

	// select RGB file and split channels for separate analysis
	run("Split Channels");

	// hop1 analysis
	close(title + " (blue)"); // close dapi channel
	close(title + " (red)"); // close zip1 channel
	selectWindow(title + " (green)"); // select hop1 channel
	saveAs("Tiff", output + File.separator + "Hop1_no_threshold" + File.separator + title0);
	setAutoThreshold("Moments dark");
		// 11/22/2022 update: I changed the Threshold to Moments. It matches Hop1 channel better
		// I am applying it now to old RGB images too in the final graphs
	setOption("BlackBackground", true); // need to clarify otherwise inverted LUT will be applied - white background
	run("Analyze Particles...", "size=0.0-Infinity display clear summarize add");
	
	// I want to omit convert to mask step before measuring foci b/c binary omits the range of signal
	// but I convert to mask after in order to save the image b/c otherwise it just saves the original nonthreshold image
	run("Convert to Mask");
	saveAs("Tiff", output + File.separator + "Hop1_threshold" + File.separator + title1);
	
	selectWindow("Results");
	saveAs("results", output + File.separator + "individual" + File.separator + title2);

	dir = getDirectory("home");
	close(title1 + ".tif");
	close(title2);
}


selectWindow("Summary");
saveAs("results",  output + File.separator + "Hop1FociArea.csv");
close("Hop1FociArea.csv");
close("Results");
close("ROI Manager");

