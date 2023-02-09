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


// Reminder: when choosing input & output, make sure file extension is .dv, not .tif like for other macro analyses scripts

function processFile(input, output, file) {
	// Do the processing here by adding your own code.
	// Leave the print statements until things work, then remove them.
	path = input + File.separator;
	full_path = path + file;
	open(full_path); // image is now open

	// set up variables
	title = getTitle();
	title1 = File.nameWithoutExtension + "_RGB";


	// if you don't select anything and press OK there will be an error
	// ** fix this w selectionType() function
	
	// if you exit out of the dialog box, there will be an error
	// waitForUser("Draw outline of yeast cell, then hit OK"); 

	// select freehand tool

	// run("Clear Outside", "stack");
	run("Split Channels");

	waitForUser("Choose best Z-stack for each channel, then hit OK");

	// duplicating best Z-stacks
	selectWindow("C1-" + title);
	run("Duplicate...", "use");
	titleBlue = getTitle();
	
	selectWindow("C2-" + title);
	run("Duplicate...", "use");
	titleGreen = getTitle();
	
	selectWindow("C3-" + title);
	run("Duplicate...", "use");
	titleRed = getTitle();

	close("C1-" + title);
	close("C2-" + title);
	close("C3-" + title);

	// In merge channels, C1 = c3/3 (red) , C2 = c2/3 (green) and C3 = c1/3 (blue dapi)
	run("Merge Channels...", "red=[" + titleRed + "] green=[" + titleGreen + "] blue=[" + titleBlue + "] create");
	
	run("Stack to RGB");
	close("Composite");
	selectWindow("Composite (RGB)");
	saveAs("Tiff", output + File.separator + title1);
	close(title1 + ".tif");
	
}
