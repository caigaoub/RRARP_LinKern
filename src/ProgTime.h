#include <chrono>

using namespace std;

class ProgTime{
public:
	chrono::system_clock::time_point   			_start;
	chrono::system_clock::time_point			_end;
	double										_elapsed_secs;
	void start_prog(){
		_start = chrono::system_clock::now();
	};

	void end_prog(){
		_end = chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = _end - _start;
		_elapsed_secs = elapsed_seconds.count();
		// return elapsed_seconds;
	};
};