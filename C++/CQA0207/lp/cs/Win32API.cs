////////////////////////////////////////////////////////////////
// MSDN Magazine -- July 2002
// If this code works, it was written by Paul DiLascia.
// If not, I don't know who wrote it.
// Compiles with Visual Studio 6 or .NET.
// Runs in Windows XP and probably Windows 2000 too.
//
using System;
using System.Text;
using System.Runtime.InteropServices;

//////////////////
// namspace to wrap Win32 API functions. Add them here as you need...
//
namespace Win32API {
	public class Win32 {
		[DllImport("user32.dll")]
		public static extern int GetClassName(int hwnd,
			StringBuilder buf, int nMaxCount);

		[DllImport("user32.dll")]
		public static extern int GetWindowText(int hwnd,
			StringBuilder buf, int nMaxCount);
	}
}
