//---------------------------------------------------------------------------
//
//    FCST: Fuel Cell Simulation Toolbox
//
//    Copyright (C) 2013 by Energy Systems Design Laboratory, University of Alberta
//
//    This software is distributed under the MIT License.
//    For more information, see the README file in /doc/LICENSE
//
//    - Description: Logging objects to be used throughout fcst
//    - Developers: Philip Wardlaw
//
//---------------------------------------------------------------------------

#include <deal.II/base/logstream.h>
#include <deal.II/base/utilities.h>

using namespace dealii;

#ifndef _FCST_LOGGING
#define _FCST_LOGGING

namespace FcstUtilities
{
    // Needed to declare extern file below
    class FCSTLogStream;    
    /**
     * Object used to output data to file and, if file attached recorded
     * to a file as well.
     */
    extern FCSTLogStream log;
    
     /**
     * The object FcstUtilities::log should be used throughout OpenFCST for console logging.
     * 
     * When OpenFCST is compiled with PETSC/MPI the logging class is changed to accommodate logging
     * from multiple processes.
     * 
     * @note This class is inspired in LogStream from deal.II
     * 
     * @note Currently, the class does not implement the push and pop functions for adding prefixes
     * 
     * @note TODO I think I could remove all pre-compiler directives by storing a pointer to a number that is initialized in the appplication
     * and that keeps track of the processor... 
     * 
     * @author M. Secanell
     */   
    class FCSTLogStream
    {
    public:
        /**
         * Standard constructor, since we intend to provide an object
         * <tt>deallog</tt> in the library. Set the standard output stream to
         * <tt>std::cerr</tt>.
         */
        FCSTLogStream (std::ostream &stream);
        
        /**
         * Destructor.
         */
        ~FCSTLogStream();
        
        /**
         * Enable output to a second stream, @p o, usually an output file.
         *
         * The optional argument @p print_job_id is not used and it is maintained here to keep
         * the same format as deal LogStream.
         * 
         * <h3> Usage </h3>
         * See SimulatorBuilder<dim>::open_logfile
         * 
         */
        void attach (std::ostream &o,
                     const bool    print_job_id = true);
        
        
        /**
         * Disable output to the second stream. You may want to call
         * <tt>close</tt> on the stream that was previously attached to this
         * object.
         */
        void detach ();        
        
        /**
         * Setup the logstream for regression test mode.
         *
         * This sets the parameters #double_threshold, #float_threshold, and
         * #offset to nonzero values. The exact values being used have been
         * determined experimentally and can be found in the source code.
         *
         * Called with an argument <tt>false</tt>, switches off test mode and
         * sets all involved parameters to zero.
         */
        /*void test_mode (bool on=true)
        {
          //  internal_log.test_mode(on);
        }
        */
        
        /**
         * Gives the default stream (<tt>std_out</tt>).
         */
        /*
        std::ostream &get_console ()
        {
           // return internal_log.get_console();
        }
        */
        
        
        /**
         * Gives the file stream.
         */
        std::ostream &get_file_stream ()
        {
           // return internal_log.get_file_stream();
        }
        
        
        /**
         * @return true, if file stream has already been attached.
         */
        bool has_file () const
        {
          return (file != 0);
        }
        
        /**
         * Push another prefix on the stack. Prefixes are automatically separated
         * by a colon and there is a double colon after the last prefix.
         *
         * A simpler way to add a prefix (without the manual need to add the
         * corresponding pop()) is to use the Prefix class.
         */
        void push (const std::string &text)
        {
            
        }
        
        
        /**
         * Remove the last prefix added with push().
         */
        void pop ()
        {
            
        }
        
        
        /**
         * Maximum number of levels to be printed on the console. This function
         * allows to restrict console output to the upmost levels of iterations.
         * Only output with less than <tt>n</tt> prefixes is printed. By calling
         * this function with <tt>n=0</tt>, no console output will be written.
         *
         * The previous value of this parameter is returned.
         */
        unsigned int depth_console (const unsigned int n)
        {
            //return internal_log.depth_console(n);
        }
        
        
        /**
         * Maximum number of levels to be written to the log file. The
         * functionality is the same as <tt>depth_console</tt>, nevertheless,
         * this function should be used with care, since it may spoile the value
         * of a log file.
         *
         * The previous value of this parameter is returned.
         */
        unsigned int depth_file (const unsigned int n)
        {}
          
        /**
         * Output a constant something through
         * this stream. This function must be @p
         * const so that member objects of this
         * type can also be used from @p const
         * member functions of the surrounding
         * class.
         */
        template <typename T>
        const FCSTLogStream &
        operator << (const T &t) const;
        
        /**
         * Treat ostream manipulators. This
         * function must be @p const so that
         * member objects of this type can also
         * be used from @p const member functions
         * of the surrounding class.
         *
         * Note that compilers want to see this
         * treated differently from the general
         * template above since functions like @p
         * std::endl are actually overloaded and
         * can't be bound directly to a template
         * type.
         */
        const FCSTLogStream &
        operator<< (std::ostream& (*p) (std::ostream &)) const;
        
        /**
         * Determine an estimate for the memory consumption (in bytes) of this
         * object. Since sometimes the size of objects can not be determined
         * exactly (for example: what is the memory consumption of an STL
         * <tt>std::map</tt> type with a certain number of elements?), this is
         * only an estimate. however often quite close to the true value.
         */
        //std::size_t memory_consumption () const;
        
    private:
        /**
         * Reference to the stream we
         * want to write to.
         */
        std::ostream&  output_stream;
        
        /**
         * Pointer to a stream, where a copy of the output is to go to. Usually,
         * this will be a file stream.
         *
         * You can set and reset this stream by the <tt>attach</tt> function.
         */
        std::ostream  *file;
    };
    
    
    // --------------------------- inline and template functions -----------    
    template <class T>
    inline const FCSTLogStream &
    FCSTLogStream::operator<< (const T &t) const
    {
      //  #ifdef OPENFCST_WITH_PETSC
        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
     //   #endif
        {
            // Print to screen
            output_stream<<t;
            
            // and to file
            if (has_file() )
                *file<<t;
        }
        
        return *this;
    }
    
    //----
    inline const FCSTLogStream & 
    FCSTLogStream::operator<< (std::ostream& (*p) (std::ostream &)) const
    {
      //  #ifdef OPENFCST_WITH_PETSC
        if (Utilities::MPI::this_mpi_process(MPI_COMM_WORLD) == 0)
     //   #endif
        {
            output_stream<<p;
            
            // and to file
            if (has_file() )
                *file<<p;
        }
        
        return *this;
    }    

}
    
#endif //_FCST_LOGGING
