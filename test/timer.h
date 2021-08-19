#include <chrono> // for std::chrono functions
#include <mpi.h>
 
class Timer final{

private:
	// Type aliases to make accessing nested type easier
	using clock_t = std::chrono::high_resolution_clock;
	using second_t = std::chrono::duration<double, std::ratio<1> >;
	
	std::chrono::time_point<clock_t> m_beg;
    

    double m_total {};
	bool running{false};
 
public:
	Timer()
        : m_beg(clock_t::now()), running{true} {}
    
    void restart(){
        m_beg = clock_t::now();
		running = true;
    }

	void pause(){
		checkpoint();
		running = false;
	}

	void stop(){
		checkpoint();
		running = false;
	}
	
	void reset(){
		m_total = 0.0;
		running = true;
		m_beg = clock_t::now();
	}

    void checkpoint(){
		if(running)
        	m_total += std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
        m_beg = clock_t::now();
    }
	
	double elapsed()
	{	
		if(running)
			stop();
		return m_total;
	}
};


class TimerMPI final{

private:
	double m_beg{};
    double m_total {};
	bool running{false};
 
public:
	TimerMPI()
        : m_beg{MPI_Wtime()}, running{true} {}
    
    void restart(){
        m_beg = MPI_Wtime();
		running = true;
    }

	void pause(){
		checkpoint();
		running = false;
	}

	void stop(){
		checkpoint();
		running = false;
	}
	
	void reset(){
		m_total = 0.0;
		running = true;
		m_beg = MPI_Wtime();
	}

    void checkpoint(){
		if(running)
        	m_total += (MPI_Wtime() - m_beg);
        m_beg = MPI_Wtime();
    }
	
	double elapsed()
	{	
		if(running)
			stop();
		return m_total;
	}
};

