////////////////////////////////////////////////////////////////
// MSDN Magazine -- July 2002
// If this code works, it was written by Paul DiLascia.
// If not, I don't know who wrote it.
// Compiles with Visual Studio 6 or .NET.
// Runs in Windows XP and probably Windows 2000 too.
//
// lp -- list processes from command line.
//
using System;
using System.Text;
using System.Diagnostics;
using System.Runtime.InteropServices;
using Win32API; // home-brew wrapper for native API

class MyApp {
	// global command-line switches
	static bool bBare = false;			// don't display header
	static bool bClassName = false;	// show window class name
	static bool bTitle = false;		// show title

	[STAThread]
	// main entry point
	static int Main(string[] args) {
		// Parse command line. Switches can come in any order.
		//
		for (int i=0, len=args.GetLength(0); i<len; i++) {
			if (args[i].StartsWith("/") || args[i].StartsWith("-") ) {
				for (int j=1; j<args[i].Length; j++) {
					switch (args[i][j]) {
					case 'b':
						bBare = true;
						break;
					case 'c':
						bClassName = true;
						break;
					case 't':
						bTitle = true;
						break;
					case '?':
					default:
						help();
						return 0;
					}
				}
			}
		}
		
		// Iterate over all processes.
		// Process.GetProcesses returns an array of them.
		//
		Int32 count=0;
		Process[] procs = Process.GetProcesses();
		bool bFirstModule=true;
		for (int i=0, len=procs.GetLength(0); i<len; i++) {
			Process p = procs[i];
			if (p.Id!=0) {
				int hwnd = p.MainWindowHandle.ToInt32();
				if (hwnd!=0) { // if has a main window:
					ProcessModule pm = p.MainModule;
					String modname = pm.ModuleName;
					if (bFirstModule) {
						if (!bBare) {
							Console.WriteLine("PID  Module       HWND     {0}{1}{2}",
								bClassName ? "class name" : "",
								bClassName && bTitle ? "/" : "",
								bTitle ? "title" : "");
							Console.WriteLine("---- ------------ -------- ----------------");
						}
						bFirstModule=false;
					}
					StringBuilder cname = new StringBuilder(256);
					Win32.GetClassName(hwnd, cname, cname.Capacity);
					StringBuilder title = new StringBuilder(256);
					Win32.GetWindowText(hwnd, title, title.Capacity);

					Console.WriteLine("{0:x4} {1}{2:x8} {3} {4}",
						p.Id, modname.PadRight(13),
						hwnd,
						bClassName ? cname.ToString() : "",
						bTitle ? title.ToString() : "");
					count++;
				}
			}
		}
		if (!bBare) {
			Console.WriteLine("----");
			Console.WriteLine("{0} proceses found.",count);
		}
		return count;
	}

	static void help() {
		Console.WriteLine("lp:       List top-level proceses.");
		Console.WriteLine("          Copyright 2002 Paul DiLascia.");
		Console.WriteLine("Format:   lp [/bct?]");
		Console.WriteLine("Options:");
		Console.WriteLine(" /b(are)  no header");
		Console.WriteLine(" /c(lass) show window class names");
		Console.WriteLine(" /t(itle) show window titles");
		Console.WriteLine("");
	}
}
