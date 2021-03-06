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

imp = IJ.openImage("/home/brian/vnand/2_2v_12min.fits");
imp.show();
imp = IJ.openImage("/home/brian/vnand/2_2v_12min_X1.fits");
imp.show();
//imp3 = IJ.openImage("/home/brian/vnand/2_2v_12min_X2.fits");
//imp3.show();
//imp4 = IJ.openImage("/home/brian/vnand/2_2v_12min_X3.fits");
//imp4.show();

// Concatenate images:

// Try following line with newlines between filenames and without; this may make sed scripts easier in bash:
IJ.run("Concatenate...", "  title=[Concatenated Stacks] image1=2_2v_12min.fits image2=2_2v_12min_X1.fits image3=[-- None --]");

var imp = IJ.getImage();

imp.show();
// Convert to 16-bit imagetype:

IJ.run(imp, "16-bit", "");

// Save image as Tiff

IJ.saveAs(imp, "Tiff", "/home/brian/vnand/2_2v_12min.tif");

// Close image
IJ.log("Video saved");
imp.close();
