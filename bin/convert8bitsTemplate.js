importClass(Packages.ij.IJ);
importClass(Packages.ij.ImagePlus);
importClass(Packages.ij.ImageStack);
importClass(Packages.ij.process.ImageProcessor);
importClass(Packages.ij.gui.Toolbar);
importClass(Packages.ij.WindowManager);
importClass(Packages.ij.plugin.frame.RoiManager);
importClass(Packages.ij.gui.GenericDialog);


// Nonfunctional approach using ImageMagick:
// convert -colorspace Gray -define dcm:display-range=reset original.FIT -auto-level -format TIFF -depth 16 output.tif


// Open images:

imp = IJ.openImage("/home/brian/vnand/sample1_23min_0v(00293).tif");
imp.show();

//var imp = IJ.getImage();

//imp.show();
// Convert to 16-bit imagetype:

IJ.run(imp, "8-bit", "");

// Save image as Tiff

IJ.saveAs(imp, "Tiff", "/home/brian/vnand/sample1_23min_0v(00293)_8bit.tif");

// Close image
IJ.log("Video saved");
imp.close();

