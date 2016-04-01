importClass(Packages.ij.IJ);
importClass(Packages.ij.ImagePlus);
importClass(Packages.ij.ImageStack);
importClass(Packages.ij.process.ImageProcessor);
importClass(Packages.ij.gui.Toolbar);
importClass(Packages.ij.WindowManager);
importClass(Packages.ij.plugin.frame.RoiManager);
importClass(Packages.ij.gui.GenericDialog);


// Following line should be built by bash script:

imp = IJ.openImage("/home/brian/working/sampleStack.tif")
IJ.log("opened image");
imp.show();

//i0 = WindowManager.getCurrentWindow();
IJ.run("StackReg", "transformation=[Rigid Body]");
IJ.saveAs("Tiff", "/home/brian/working/sampleStack_bgs.tif");

IJ.log("Video saved");
imp.close();



