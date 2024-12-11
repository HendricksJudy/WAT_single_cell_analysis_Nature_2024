// Rename the current image to "image"
rename("image");

// Split the image into individual channels (e.g., RGB or multi-channel image)
run("Split Channels");

// Select the channel 2 image and close it
selectImage("C2-image");
close();

// Select the channel 1 image
selectImage("C1-image");

// Duplicate the channel 1 image for processing
run("Duplicate...", " ");

// Work on the duplicate of the channel 1 image
selectImage("C1-image-1");

// Apply sharpening filter to enhance edges
run("Sharpen");

// Apply Gaussian blur with a sigma of 5 to smooth the image
run("Gaussian Blur...", "sigma=5");

// Subtract the blurred image from the original channel 1 image to enhance structures
imageCalculator("Subtract create", "C1-image","C1-image-1");

// Select the result of the subtraction
selectImage("Result of C1-image");

// Apply sharpening to further enhance edges
run("Sharpen");

// Apply Otsu thresholding to binarize the image (dark background assumed)
setAutoThreshold("Otsu dark no-reset");
setOption("BlackBackground", true);

// Convert the binary image into a mask
run("Convert to Mask");

// Apply the watershed algorithm twice to separate touching objects
run("Watershed");
run("Watershed");

// Close the duplicate of the channel 1 image
selectImage("C1-image-1");
close();

// Rename the processed binary image to "dapi_binary"
selectImage("Result of C1-image");
rename("dapi_binary");

// Select the channel 3 image for processing
selectImage("C3-image");

// Duplicate the channel 3 image for processing
run("Duplicate...", " ");

// Apply sharpening filter to enhance edges
run("Sharpen");

// Apply Gaussian blur with a sigma of 5 to smooth the image
run("Gaussian Blur...", "sigma=5");

// Subtract the blurred image from the original channel 3 image to enhance structures
imageCalculator("Subtract create", "C3-image","C3-image-1");

// Select the result of the subtraction
selectImage("Result of C3-image");

// Apply Li thresholding to binarize the image (dark background assumed)
setAutoThreshold("Li dark no-reset");
setOption("BlackBackground", true);

// Convert the binary image into a mask
run("Convert to Mask");

// Remove small specks or noise from the binary mask
run("Despeckle");

// Apply the watershed algorithm twice to separate touching objects
// Uncomment the following line if you want additional watershed processing
// run("Watershed");
run("Watershed");

// Close the duplicate of the channel 3 image
selectImage("C3-image-1");
close();

// Rename the processed binary image to "wga_binary"
selectImage("Result of C3-image");
rename("wga_binary");

// Duplicate the "wga_binary" image for skeletonization
selectImage("wga_binary");
run("Duplicate...", " ");

// Apply skeletonization to the duplicate binary image
run("Skeletonize");

// Rename the skeletonized image to "wga_binary_skeleton"
selectImage("wga_binary-1");
rename("wga_binary_skeleton");

// Close the channel 3 image
selectImage("C3-image");
close();

// Close the channel 1 image
selectImage("C1-image");
close();