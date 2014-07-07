#include <iostream>
#include <sstream>
#include <omp.h>
//#include <unistd.h>

 
class ParallelStream{
    std::ostringstream stdStream;
public:
    ParallelStream(){}
    template <class T>
    ParallelStream& operator<<(const T& inData){
        stdStream << inData;
        return *this;
    }
    std::string toString() const{
        return stdStream.str();
    }
};


int main(){
std::cout << "Basic std::cout" << std::endl;
    //#pragma omp parallel num_threads(4)
    //{
    
    std::ostream *output = &std::cout;
    
    #pragma omp parallel for schedule(dynamic) 
    for (size_t rep_i=0; rep_i < 100; ++rep_i) {
        //std::cout << (ParallelStream() <<"I am " << omp_get_thread_num() << " working at section 3n").toString()<< std::endl;
        //std::cout << (ParallelStream() <<"I am " << omp_get_thread_num() << " working at section 3n"<<"\n").toString();// std::endl;
        ParallelStream mytmp_strm;
        mytmp_strm << "I am " << omp_get_thread_num() << " working at section 3n"<<"\n";
        mytmp_strm << "This is the second line  " << omp_get_thread_num() << " working at section 3n"<<"\n";
        *output << mytmp_strm.toString();// std::endl;

    }
 
    *output << "With Cout_MP" << std::endl;
    //#pragma omp parallel num_threads(4)
    //{
        //Cout_MP::Print() << "I want " << "to print many " << "things " << "thread : " << omp_get_thread_num() << Cout_MP::ENDL;
    //}
 
    return 0;
}
