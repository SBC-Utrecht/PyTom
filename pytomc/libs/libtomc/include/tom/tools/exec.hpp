/****************************************************************************//**
 * \file exec.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    18.02.2007
 * \brief A class to execute a child process (using fork/exec) and retrieve the standard output and error
 *
 * All functions are inline.
 * Using this file "should" be thread save.
 *******************************************************************************/
#ifndef ___INCLUDE__TOM__TOOLS__EXEC_HPP__
#define ___INCLUDE__TOM__TOOLS__EXEC_HPP__



#include <errno.h>
#include <vector>
#include <stdexcept>
#include <unistd.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/wait.h>


#include <tom/tools/snippets.hpp>




namespace tom {
namespace tools {



template<typename T> void exec(T executor, std::vector<char> &stdout, std::vector<char> &stderr, int *status);







/****************************************************************************//**
 * \brief class for handling a pipe.
 *
 * The constructor opens a pipe (see the manual page for \c pipe) and the
 * destructor closes it.
 *******************************************************************************/
class Pipe {
public:
    enum direction { READ=0, WRITE=1 };

    Pipe() {
        if (::pipe(fd) != 0) {
            const int errno_ = errno;
            throw std::runtime_error("error creating pipe: " + tom::tools::strerror(errno_));
        }
        fd_open[0] = fd_open[1] = true;
    }
    ~Pipe() {
        // The destructor should never throw an exception. Silently ignore an error.
        close_nothrow();
    }
    void close() {
        close_fd(0);
        close_fd(1);
    }
    void close_nothrow() {
        close_nothrow_fd(0);
        close_nothrow_fd(1);
    }
    void close_fd(int i) {
        // i should be either 0 or 1 (otherwise it defaults to 1).
        if (fd_open[i=(i?1:0)]) {
            fd_open[i] = false;
            if (::close(fd[i]) != 0) {
                const int errno_ = errno;
                throw std::runtime_error("error closing pipe: " + tom::tools::strerror(errno_));
            }
        }
    }
    void close_nothrow_fd(int i) {
        // i should be either 0 or 1 (otherwise it defaults to 1).
        if (fd_open[i=(i?1:0)]) {
            fd_open[i] = false;
            ::close(fd[i]);
        }
    }
    int operator[](int i) const {
        i = i ? 1 : 0;
        if (!fd_open[i]) { throw std::runtime_error("file descriptor already closed."); }
        return fd[i];
    }
    bool is_open_fd(int i) const {
        return fd_open[i ? 1 : 0];
    }
private:
    int fd[2];
    bool fd_open[2];
}; // class tom::tools::Pipe



/****************************************************************************//**
 * \brief Caller object for the template function tom::tools::exec.
 *
 * The template function tom::tools:exec needs an object with the operator ()
 * defined to perform the actuall exec system call.
 * This class provides the functionality of execvp and execv.
 *******************************************************************************/
class Execv {
public:
    Execv(const std::string &file, const std::vector<std::string> &argv, bool use_path);
    inline int operator()() const {
        if (use_path_) {
            ::execvp(&file_[0], &argv_[0]);
        } else {
            ::execv(&file_[0], &argv_[0]);
        }
        return errno;
    }
private:
    // Hide the copy constructor and the assignment operator.
    Execv(const Execv &c) { assert(0); }
    Execv &operator=(const Execv &c) { assert(0); return *this; }

    bool use_path_;
    std::vector<char> file_; // use a character vector instead of std::string because execv wants non-const pointers.
    std::vector<std::vector<char> > argv_container_; // contains the character arrays
    std::vector<char *> argv_;
}; // class tom::tools::Execv


} // namespace tools
} // namespace tom






// inline functions.



/****************************************************************************//**
 * \brief Constructor to remember the parameters for calling exec.
 *
 * \param[in] file Filename of the executable to exec.
 * \param[in] argv The input arguments of the program.
 * \param[in] use_path Decide whether to use execvp (use_path==true) or execv
 *      to search for the filename in the $PATH variable (See manual pages).
 *
 * The class keeps copies of the parameters \c file and \c argv, which can therefore
 * be altered savely after the call of the constructor.
 *******************************************************************************/
inline tom::tools::Execv::Execv(const std::string &file, const std::vector<std::string> &argv, bool use_path)
    : use_path_(use_path),
      file_(),
      argv_container_(),
      argv_() {

    if (file.empty()) {
        throw std::invalid_argument("the file parameter can not be empty.");
    }

    const std::size_t argc = argv.size();

    // Reserve the space for the parameters.
    argv_container_.resize(argc);
    argv_.resize(1 + argc);

    // Copy the filename to the internal filename buffer.
    file_.reserve(file.size()+1);
    file_.assign(file.begin(), file.end());
    file_.push_back(0);
    //argv_[0] = &file_[0]; // argv[0] points to the filename...

    // Copy the parameters.
    for (std::size_t i=0; i<argc; i++) {
        std::vector<char> &c = argv_container_[i];
        const std::string &s = argv[i];
        c.reserve(s.size()+1);
        c.assign(s.begin(), s.end());
        c.push_back(0);
        argv_[i] = &c[0];
    }

    // Set final NULL.
    argv_[argc] = 0;
}


/****************************************************************************//**
 * \brief runs an other program by calling exec and fork.
 *
 * \param[in] executor Template argument. Must be callable as function without parameters
 *    and return an error status in case of failure (preferably \c errno).
 *    If no error occures in \a executor, the function should not return as the process
 *    is replaced by a new process. Calling <c>executor()<c>  should perform the
 *    exec system call...
 * \param[out] stdout: Vector containing the standard output of the called program.
 * \param[out] stderr: Vector containing the standard error of the called program.
 * \param[out] status: Gives the return status of the program as returned by waidpid.
 *    If NULL, the status is not returned (but the process still waits for its child).
 *
 * The function forks the process and the child process calles \c exec (by invoking
 * <c>executor()<c>). In case of an error a std::runtime_error is thrown.
 * Cases of error may be for exeampl a fail to fork, an error for the child to call exec.
 *
 * \todo There is still a problem that the child could fill up the buffers for
 * stdout/stderr in the parent. The buffer constantly grows as the child outputs
 * something.
 *******************************************************************************/
template<typename T>
inline void tom::tools::exec(T executor, std::vector<char> &stdout, std::vector<char> &stderr, int *status) {

    int status_tmp;
    int &status_real = status ? *status : status_tmp;
    status_real = 0;

    stdout.resize(0);
    stderr.resize(0);
    int errno_;
    int ret;

    // Pipes to communicate with the child.
    // both for stdout, stderr and a channel to transmit the error if exec failed.
    tom::tools::Pipe fd_error, fd_stdout, fd_stderr;

    std::size_t stdout_size = 0;
    std::size_t stderr_size = 0;

    pid_t pid = fork();

    std::string error_msg = "unknown error.";

    if (pid < 0) {
        // Fork failed.
        errno_ = errno;
        fd_error.close_nothrow();
        fd_stdout.close_nothrow();
        fd_stderr.close_nothrow();
        throw std::runtime_error("could not fork: " + tom::tools::strerror(errno_));
    } else if (0 == pid) {
        // child process.
        try {

            // set the close on exec flag for the error pipe so that the parent process
            // sees when the child sucessfully called exec.
            {
                // get the flags of the file descriptor.
                int fd_ = fd_error[Pipe::WRITE];
                if ((ret = fcntl(fd_, F_GETFD)) == -1) {
                    errno_ = errno;
                    throw std::runtime_error("Error calling fcntl(fd, F_GETFD) on the error pipe: " + tom::tools::strerror(errno_));
                }

                // make the write end close-on-exec
                int fd_flags = ret | FD_CLOEXEC;
                if (fcntl(fd_, F_SETFD, fd_flags) == -1) {
                    errno_ = errno;
                    throw std::runtime_error("Error calling fcntl(fd, F_SETFD, flags) on the error pipe: " + tom::tools::strerror(errno_));
                }
            }
            fd_error.close_fd(Pipe::READ);

            // Close the STDIN of the cild.
            if (close(STDIN_FILENO) != 0) {
                errno_ = errno;
                throw std::runtime_error("Could not close STDIN: " + tom::tools::strerror(errno_));
            }

            // connect the STDOUT of the child with the parent.
            if (dup2(fd_stdout[Pipe::WRITE], STDOUT_FILENO) == -1) {
                errno_ = errno;
                throw std::runtime_error("Error redirect STDOUT to parent: " + tom::tools::strerror(errno_));
            }
            fd_stdout.close();

            // connect the STDERR of the child with the parent.
            if (dup2(fd_stderr[Pipe::WRITE], STDERR_FILENO) == -1) {
                errno_ = errno;
                throw std::runtime_error("Error redirect STDERR to parent: " + tom::tools::strerror(errno_));
            }
            fd_stderr.close();

            errno = 0;
            try {
                errno_ = executor();
            } catch (std::exception &e) {
                errno_ = errno;
                throw std::runtime_error(std::string("Error calling exec: \"") + e.what() + "\".");
            } catch (...) {
                errno_ = errno;
                throw std::runtime_error("Error calling exec (unknown exception): \"" + tom::tools::strerror(errno_) + "\".");
            }
            error_msg = "exec failed (" + tom::tools::strerror(errno_) + ")";
        } catch (std::exception &e) {
            error_msg = e.what();
        } catch (...) {
            error_msg = "Unexpected exception in child process.";
        }
        error_msg = "child process: " + error_msg;
        // try to write the error message to the parent.
        try {
            write(fd_error[Pipe::WRITE], error_msg.c_str(), error_msg.size()+1);
        } catch (...) {
        }
        // che child exits with error code -1.
        exit(-1);
    }


    try {
        // parent process.
        pid_t wpid;

        const std::size_t BUFFSIZE = 1024;
        ssize_t cread = 0;
        try {
            // Close the write end of the pipes.
            fd_error.close_fd(Pipe::WRITE);
            fd_stdout.close_fd(Pipe::WRITE);
            fd_stderr.close_fd(Pipe::WRITE);

            {
                // Read the error message from the error-pipe (where the child tries
                // to send its error messages that happend before or while calling exec.
                int fd = fd_error[Pipe::READ];
                std::vector<char> cbuff(BUFFSIZE);
                std::size_t cbuff_free = BUFFSIZE;
                cread = 0;
                do {
                    cbuff_free -= cread;
                    if (!cbuff_free) {
                        cbuff.resize(cbuff.size() + BUFFSIZE, 0);
                        cbuff_free = BUFFSIZE;
                    }
                    do {
                        cread = read(fd, &cbuff[cbuff.size() - cbuff_free], cbuff_free);
                    } while (cread==-1 && errno==EINTR); // In case of returned by a signal: repeat reading.
                } while (cread > 0);
                if (cread == -1) {
                    errno_ = errno;
                    throw std::runtime_error(std::string("error reading the error message fron the child: ") + tom::tools::strerror(errno_));
                }
                if (cbuff.size() > cbuff_free) {
                    // The buffer is already null-terminated because resize sets the content.
                    // Morover there are still cbuff_free (>0) elements which do not belong to the message (and are therefore 0).
                    throw std::runtime_error(&cbuff[0]);
                }
                fd_error.close();
            }

            // At this point, the child process successfully called exec.
            // Now read the STDOUT and STDERR...
            fd_set rfds;
            bool fd_stdout_open = true;
            bool fd_stderr_open = true;
            int fd_stdout_errno = 0;
            int fd_stderr_errno = 0;
            int fd_stderr_ = fd_stderr[Pipe::READ];
            int fd_stdout_ = fd_stdout[Pipe::READ];
            int nfds = (fd_stderr_>fd_stdout_ ? fd_stderr_ : fd_stdout_) + 1; // Get the number of the highest file decriptor for select.
            int i;
            // Now read from the pipes into the vectors.
            // Use select to read only when there is really data available.
            // Otherwise it likely will deadlock.

            // Repeat as long there are open streams (pipes).
            do {
                // Init the bit-mask.
                FD_ZERO(&rfds);
                if (fd_stdout_open) { FD_SET(fd_stdout_, &rfds); }
                if (fd_stderr_open) { FD_SET(fd_stderr_, &rfds); }
                if ((i = select(nfds, &rfds, NULL, NULL, NULL)) < 0) {
                    // error in select.
                    errno_ = errno;
                    if (i==-1 && errno_ == EINTR) { continue; }
                    throw std::runtime_error(std::string("Error while waiting to read from streams (select): ") + tom::tools::strerror(errno_));
                }
                if (FD_ISSET(fd_stdout_, &rfds)) {
                    // Read from stream.
                    stdout.resize(stdout_size+BUFFSIZE, 0);
                    cread = read(fd_stdout_, &stdout[stdout_size], BUFFSIZE);
                    if (cread == -1) {
                        errno_ = errno;
                        if (errno != EINTR) {
                            // silently ignore this error for the moment. Throw an exception later by checking fd_stdout_errno
                            fd_stdout_open = false;
                            fd_stdout_errno = errno_;
                            fd_stdout.close_nothrow();
                        }
                    } else if (cread == 0) {
                        // EOF
                        fd_stdout_open = false;
                        fd_stdout.close_nothrow();
                    } else {
                        stdout_size += cread;
                    }
                }
                if (FD_ISSET(fd_stderr_, &rfds)) {
                    // Read from stream.
                    stderr.resize(stderr_size+BUFFSIZE, 0);
                    cread = read(fd_stderr_, &stderr[stderr_size], BUFFSIZE);
                    if (cread == -1) {
                        errno_ = errno;
                        if (errno_ != EINTR) {
                            // silently ignore this error for the moment. Throw an exception later by checking fd_stderr_errno
                            fd_stderr_open = false;
                            fd_stderr_errno = errno_;
                            fd_stderr.close_nothrow();
                        }
                    } else if (cread == 0) {
                        // EOF
                        fd_stderr_open = false;
                        fd_stderr.close_nothrow();
                    } else {
                        stderr_size += cread;
                    }
                }
            } while (fd_stdout_open || fd_stderr_open);

            if (fd_stdout_errno!=0 && fd_stderr_errno!=0) {
                throw std::runtime_error("Error reading STDOUT (\"" + tom::tools::strerror(fd_stdout_errno) + "\") and STDERR (\"" + tom::tools::strerror(fd_stderr_errno) + "\") from child process.");
            } else if (fd_stdout_errno != 0) {
                throw std::runtime_error("Error reading STDOUT from child process: \"" + tom::tools::strerror(fd_stdout_errno) + "\".");
            } else if (fd_stderr_errno != 0) {
                throw std::runtime_error("Error reading STDERR from child process: \"" + tom::tools::strerror(fd_stderr_errno) + "\".");
            }
        } catch (...) {
            // In any case wait for the termination of the child process...
            do {
                wpid = waitpid(pid, &status_real, 0);
            } while (wpid==-1 && errno==EINTR); // if a signal woke the parent, sleep again....

            // Rethorw the original exception.
            throw;
        }

        // wait for child process.
        do {
            wpid = waitpid(pid, &status_real, 0);
        } while (wpid==-1 && errno==EINTR); // if a signal woke the parent, sleep again....
        if (wpid == -1) {
            errno_ = errno;
            throw std::runtime_error("Error waiting for child (waitpid) \"" + tom::tools::strerror(errno) + "\".");
        }
    } catch (...) {
        stdout.resize(stdout_size);
        stderr.resize(stderr_size);
        throw;
    }
    stdout.resize(stdout_size);
    stderr.resize(stderr_size);
}




#endif










