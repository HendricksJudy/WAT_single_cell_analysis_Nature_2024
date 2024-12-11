// Run "Analyze Particles" before this script. Use on image obtained from "xenium_get_wga_binary_final.ijm"
// Suggested parameters for Analyze Particles: 
// - Minimum size: 5 for nuclei, 30 for wga, 0 for wga_skeleton. 
// WGA was used for Area correlation analysis and transcript assignment.
// After running analyze particles, run Measure to get descriptors such as Area, as well as the centroid position for each object

finalx = newArray(); // Initialize an array to store the x-coordinates of points.
finaly = newArray(); // Initialize an array to store the y-coordinates of points.
finalCell = newArray(); // Initialize an array to store the ROI (Region of Interest) indices.

numROIs = roiManager("count"); // Get the total number of ROIs in the ROI Manager.

for (i = 0; i < numROIs; i++) { 
    roiManager("Select", i); // Select the ROI at index `i` from the ROI Manager.
    
    Roi.getContainedPoints(xpoints, ypoints); // Get the x and y coordinates of all points contained within the selected ROI.
    
    name = Array.copy(xpoints); // Copy the x-coordinate array to use as a template for ROI indexing.
    name = Array.fill(name, i + 1); // Replace all elements in the copied array with the ROI index (i + 1).
    
    finalx = Array.concat(finalx, xpoints); // Concatenate the x-coordinates of the current ROI with the overall x-coordinate array.
    finaly = Array.concat(finaly, ypoints); // Concatenate the y-coordinates of the current ROI with the overall y-coordinate array.
    finalCell = Array.concat(finalCell, name); // Concatenate the ROI indices for the current points with the overall indices array.
    
    run("Collect Garbage"); // Free up memory by running garbage collection.
}

Table.create("_ij_result"); // Create a new Results Table named "_ij_result".

Table.setColumn("object", finalCell); // Add the ROI indices as a column named "object".
Table.setColumn("xpoints", finalx); // Add the x-coordinates as a column named "xpoints".
Table.setColumn("ypoints", finaly); // Add the y-coordinates as a column named "ypoints".

run("Collect Garbage"); // Run garbage collection again to free memory.