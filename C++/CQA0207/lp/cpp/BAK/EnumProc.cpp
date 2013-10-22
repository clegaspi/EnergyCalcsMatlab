////////////////////////////////////////////////////////////////
// MSDN Magazine -- July 2002
// If this code works, it was written by Paul DiLascia.
// If not, I don't know who wrote it.
// Compiles with Visual Studio 6.0
// Runs in Windows XP and probably Windows 2000 too.
//
#include "stdafx.h"
#include "EnumProc.h"

CProcessIterator::CProcessIterator()
{
	m_pids = NULL;
}

CProcessIterator::~CProcessIterator()
{
	delete [] m_pids;
}

//////////////////
// Get first process.
// Call EnumProcesses to initialize entire array, then return first one.
//
DWORD CProcessIterator::First()
{
	m_current = (DWORD)-1;
	m_count = 0;
	DWORD nalloc = 1024;
	do {
		delete [] m_pids;
		m_pids = new DWORD [nalloc];
		if (EnumProcesses(m_pids, nalloc*sizeof(DWORD), &m_count)) {
			m_count /= sizeof(DWORD);
			m_current = 1;						 // skip IDLE process
		}
	} while (nalloc <= m_count);

	return Next();
}

CProcessModuleIterator::CProcessModuleIterator(DWORD pid)
{
	m_hProcess = OpenProcess(PROCESS_QUERY_INFORMATION | PROCESS_VM_READ,
		FALSE, pid);
}

CProcessModuleIterator::~CProcessModuleIterator()
{
	CloseHandle(m_hProcess);
	delete [] m_hModules;
}

//////////////////
// Get first module in process. Call EnumProcesseModules to
// initialize entire array, then return first module.
//
HMODULE CProcessModuleIterator::First()
{
	m_count = 0;
	m_current = (DWORD)-1; 
	m_hModules = NULL;
	if (m_hProcess) {
		DWORD nalloc = 1024;
		do {
			delete [] m_hModules;
			m_hModules = new HMODULE [nalloc];
			if (EnumProcessModules(m_hProcess, m_hModules,
				nalloc*sizeof(DWORD), &m_count)) {
				m_count /= sizeof(HMODULE);
				m_current = 0;
			}
		} while (nalloc <= m_count);
	}
	return Next();
}

//////////////////
// Generic top-level window iterator encapsulates ::EnumWindows
//
CWindowIterator::CWindowIterator(DWORD nAlloc)
	: m_current(0), m_count(0)
{
	assert(nAlloc>0);
	m_hwnds = new HWND [nAlloc];
	m_nAlloc = nAlloc;
}

CWindowIterator::~CWindowIterator()
{
	delete [] m_hwnds;
}

//////////////////
// Enumerator fn: pass to virtual fn
//
BOOL CALLBACK CWindowIterator::EnumProc(HWND hwnd, LPARAM lp)
{
	return ((CWindowIterator*)lp)->OnEnumProc(hwnd);
}

//////////////////
// Virtual enumerator proc: add HWND to array if OnWindows is TRUE.
//
BOOL CWindowIterator::OnEnumProc(HWND hwnd)
{
	if (OnWindow(hwnd)) {
		if (m_count < m_nAlloc)
			m_hwnds[m_count++] = hwnd;
	}
	return TRUE; // keep looking
}

//////////////////
// Process window iterator: special case to iterate main windows of a process
//
CMainWindowIterator::CMainWindowIterator(DWORD pid, BOOL bVis,
	DWORD nAlloc) : CWindowIterator(nAlloc), m_pid(pid), m_bVisible(bVis)
{

}

CMainWindowIterator::~CMainWindowIterator()
{

}

//////////////////
// virtual override: is this window a main window of my process?
//
BOOL CMainWindowIterator::OnWindow(HWND hwnd)
{
	if (!m_bVisible || (GetWindowLong(hwnd,GWL_STYLE) & WS_VISIBLE)) {
		DWORD pidwin;
		GetWindowThreadProcessId(hwnd, &pidwin);
		if (pidwin==m_pid)
			return TRUE;
	}
	return FALSE;
}
