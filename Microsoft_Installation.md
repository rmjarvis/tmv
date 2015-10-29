﻿#summary Instructions for building TMV in Microsoft Visual C++

# Installing using Microsoft Visual C++ #

Using TMV with Microsoft Visual C++ is a bit different from the other compilers, since it has a Windows framework for building programs, rather than a command line.

(There is a way to do compile on the command line, but I suspect that will not be the usual way that people will want to use the TMV library.  If you are interested in compiling this on the command line, I did have success installing SCons and then, within the Visual Studio 2008 Command Prompt, using the command: `C:\Python26\Scripts\scons CXX=cl IMPORT_ENV=true`.)

1. Download the tarball as described on the [Installation](Installation.md) page.

There are many Windows utilities that can unpack the tarball.  With [IZArc](http://www.izarc.org/), for example, you right click on the file `tmv0.70.tar.gz` in Windows Explorer, select IZArc, then select Extract Here.  This should make a directory called `tmv0.70` which has all of the source code and other files for the TMV library.

2. Start Microsoft Visual C++.

I have Visual C++ 2008 Express Edition, so all instructions below about menus and such refer that that edition.  I would assume that other editions have their menus arranged similarly, but there may be some differences.

3. Open the "Solution" file `tmvtest1.sln`.

Each test suite needs to be made individually, so I'll give detailed instructions for the first test program.  The other two are made the same way.

Go to the File menu.  Select Open.  Then select Project/Solution...  Look for the file `tmvtest1.sln` and select it.

This includes a project for the first test program, `tmvtest1.vcproj`, and a project for the main TMV library, `tmv.vcproj`.

4. Select Debug or Release mode.

There are two modes for building programs in Visual C++: Debug and Release.   You can choose either one.  Dubug mode will compile faster and  Release mode will execute faster.

When writing your own code, you will probably want to start with Debug mode and then switch to Release when everything is working, so you will end up compiling the TMV library in both modes anyway.

You can select the mode in a pull down menu in the top button bar.  With my setup (which I think is the default), it is directly to the right of a green arrow.

5. Build `tmvtest1`.

Go to the Build menu.  Select Build Solution.

6. Run `tmvtest1`.

Go to the Debug menu.  Select Start Without Debugging.  (Of course, you can instead choose Start Debugging if you'd rather, but that will make it run slower.)

A console window should open up to show the output of the program.  When it is done, you should see the message,
"Press any key to continue . . . "

The previous output lines should all read _Something_` passed all tests`. If the last line starts with `Error`,  then please post an issue at http://code.google.com/p/tmv-cpp/issues/list about the problem.

7. Repeat for `tmvtest2` and `tmvtest3`.

Note: the solution file for `tmvtest2` includes the project for the TMV library with symmetric and banded matrices, `tmv_symband.vcproj`.

8. Include `tmv.vcproj` (and `tmv_symband.vcproj` if necessary) in the solution for your own project.

Any program that uses the TMV library needs to include the project `tmv.vcproj`.  To do this, Go to the File, select Add, Existing Project... In the `tmv\tmvversion` directory, look for the file `tmv.vcproj` and select it.  If you are going to be using symmetric and/or banded matrices, then you should also select `tmv_symaband.vcproj`.

Then, in the Solution Explorer window, select your project.  Then go to the Project menu and select Properties.  Select Common Properties, then Framework and References.   Click the Add New Reference... button.  Select the TMV library (or libraries), and press OK.

Next you need to tell your program where to look for the TMV header files.  So on the Property Page, select All Configurations at the top.  Then go to Configuration Properties, C/C++, General. At the entry for Additional Include Directories, click the `[...]` button at the right.  Click on the folder icon to add a new line, and then click the `[...]` button.  Browse to the  `tmv\tmvversion` directory, and select the include directory.  Press OK three times.

Now you should be all set to use the TMV library in your program.