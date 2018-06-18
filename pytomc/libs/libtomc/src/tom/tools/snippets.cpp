/***********************************************************************//**
 * \file snippets.cpp
 * \brief Implementations of snippets.hpp
 * \author  Thomas Haller
 * \version 0.2
 * \date    05.02.2007
 *
 * Implementation of snippets.hpp
 **************************************************************************/
#include <tom/tools/snippets.hpp>


#include <set>
#include <map>
#include <ostream>
#include <iomanip>

#include <cstring>
#include <cerrno>
#include <ctime>
#include <cmath>

#include <unistd.h>

#include <boost/lexical_cast.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/format.hpp>




/****************************************************************************//**
 * Uses POSIX nanosleep to sleep for a (non integer) number of seconds.
 *
 * Does not check for interrupted sleeps etc. This is only for convenience
 * to wrapp setting the struct timespec.
 *******************************************************************************/
void tom::tools::nanosleep(double sec) {
    if (sec > 0) {
        struct timespec t;
        t.tv_sec = static_cast<time_t>(trunc(sec));
        t.tv_nsec = static_cast<long>(1000000000. * (sec - trunc(sec)));
        nanosleep(&t, 0);
    }
}





/****************************************************************************//**
 * Helper function that removes spaces at the start and end of a string.
 *******************************************************************************/
void tom::tools::str_trim_spaces(std::string &str) {

    // Trim Both leading and trailing spaces
    size_t startpos = str.find_first_not_of(" \t"); // Find the first character position after excluding leading blank spaces
    size_t endpos = str.find_last_not_of(" \t"); // Find the first character position from reverse af

    // if all spaces or empty return an empty string
    if((std::string::npos == startpos) || (std::string::npos == endpos)) {
        str = "";
    } else {
        str = str.substr(startpos, endpos-startpos+1);
    }
}



/****************************************************************************//**
 * \brief Returns the current working directory.
 *******************************************************************************/
std::string tom::tools::getcwd() {
    std::vector<char> buf(128, 0);
    char *pbuf;
    int errno_;
    do {
        pbuf = ::getcwd(&buf[0], buf.size());
        errno_ = errno;
        if (!pbuf && errno_==ERANGE) {
            // The buffer is too small...
            buf.resize(buf.size()*2, 0);
        }
    } while (!pbuf && errno_==ERANGE);

    return &buf[0];
}

/****************************************************************************//**
 * \brief Returns the current time as a string
 *
 * The time format is chosen in order to allow an easy pasing of the time.
 *******************************************************************************/
std::string tom::tools::now2str(const std::string &s) {
    time_t now = time(NULL);
    std::ostringstream ss;
    struct tm now_tm;
    localtime_r(&now, &now_tm);
    if (s == "only_numbers") {
        ss << (now_tm.tm_year+1900) <<
            std::setfill ('0') << std::setw(2) << (now_tm.tm_mon+1) <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_mday << '_' <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_hour <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_min <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_sec;
    } else if (s == "YYYYMMDD-HHMMSS") {
        ss << (now_tm.tm_year+1900) <<
            std::setfill ('0') << std::setw(2) << (now_tm.tm_mon+1) <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_mday << '-' <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_hour <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_min <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_sec;
    } else {
        ss << (now_tm.tm_year+1900) << "." <<
            std::setfill ('0') << std::setw(2) << (now_tm.tm_mon+1) << "." <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_mday << "," <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_hour << ":" <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_min << ":" <<
            std::setfill ('0') << std::setw(2) << now_tm.tm_sec;
    }
    return ss.str();
}



namespace {
template<typename T> char *tom__tools__strerror(T v, char *buf);
template<> inline char *tom__tools__strerror(char *v, char *buf) { return v  ; }
template<> inline char *tom__tools__strerror(int   v, char *buf) { return buf; }
}
/****************************************************************************//**
 * \brief Wrappes strerror_r.
 *******************************************************************************/
std::string tom::tools::strerror(int errno_) {
    #define LENGTH 2048
    char buf[LENGTH];
    char *p = ::tom__tools__strerror(::strerror_r(errno_, buf, LENGTH), buf);
    buf[LENGTH-1] = 0;
    return p;
}


/****************************************************************************//**
 * \brief Insert a prefix for each line in the string.
 *******************************************************************************/
std::string tom::tools::prefix_lines(const std::string &s, const std::string &prefix) {

    if (prefix.empty()) {
        return s;
    }
    if (s.empty()) {
        return prefix;
    }

    std::ostringstream ss;
    std::vector<char> buf;
    const std::size_t ssize = s.size();
    buf.reserve(ssize+1);
    buf.insert(buf.begin(), s.begin(), s.end());
    buf.push_back(0);
    char *pbuf = &buf[0];
    char *pend = pbuf + (ssize+1);
    char *idx;
    while (pbuf < pend) {
        if ((idx = std::strchr(pbuf, '\n'))) {
            idx[0] = 0;
            ss << prefix << pbuf << '\n';
            pbuf = idx+1;
        } else {
            if (pbuf[0]) {
                ss << prefix << pbuf;
            }
            break;
        }
    }
    return ss.str();
}



/****************************************************************************//**
 *
 *******************************************************************************/
std::string tom::tools::split_line(std::string &s, std::size_t max_line_length) {
    #define DELIMITER " \n\t"
    std::string s0;
    if (max_line_length < 1) {

    } else if (max_line_length >= s.size()) {
        s0.swap(s);
        s.clear();
    } else {
        std::string s_;
        s.substr(0, max_line_length+1).swap(s_);
        std::string::size_type idx = s_.find_last_of(DELIMITER);
        if (idx == std::string::npos) {
            idx = s.find_first_of(DELIMITER, max_line_length-1);
        }

        if (idx != std::string::npos) {
            s.substr(0, idx).swap(s0);
            s.erase(0, idx+1);
        } else {
            s0.swap(s);
            s.clear();
        }
    }

    return s0;
    #undef DELIMITER
}

/****************************************************************************//**
 *
 *******************************************************************************/
void tom::tools::print_parameters(std::ostream &clog, const std::string &name, const std::string &text, const std::string &prefix, std::size_t max_line_length) {

    std::string s(text);
    std::size_t i;
    if (!name.empty()) {
        i = name.size() + 2;
        clog << name << ": " << tom::tools::split_line(s, i>max_line_length ? 0 : max_line_length-i) << "\n";
    }

    i = prefix.size();
    if (max_line_length < i) {
        max_line_length = 1;
    } else {
        max_line_length -= i;
    }

    while (!s.empty()) {
        clog << prefix << tom::tools::split_line(s, max_line_length) << "\n";
    }

}




/****************************************************************************//**
 * \param[in] instr The input string to be parsed to an integer range.
 *   The single elements are separated by ','. '-' specifies a range (with both
 *   numbers included. Negative numbers are not allowed. Spaces are not allowed.
 * \param[out] range The found numbers. Sorted in accending order.
 * \return True if the string could be successfully parsed. Otherwise false.
 *
 * Each value comes only one time in the result. Empty ranges occure never.
 * The output \a range will be sorted.
 *******************************************************************************/
template<typename T>
bool tom::tools::str2range(const std::string &instr, std::vector<T> &range, char delimiter, char delimiter_range) {

    range.resize(0);
    if (instr.empty()) {
        return false;
    }

    std::set<T> srange;
    std::string temp;
    std::string str = instr;

    std::string::size_type pos;
    T val1, val2, i;
    try {
        do {
            if ((pos=str.find(delimiter, 0)) == std::string::npos) {
                pos = str.size();
            }
            temp = str.substr(0, pos);
            str.erase(0, pos + 1);
            pos = temp.find(delimiter_range, 0);
            if (pos == std::string::npos) {
                srange.insert(boost::lexical_cast<T>(temp));
            } else if (temp.find(delimiter_range, pos+1) == std::string::npos) {
                val1 = boost::lexical_cast<T>(temp.substr(0, pos));
                temp.erase(0, pos+1);
                val2 = boost::lexical_cast<T>(temp);
                for (i=val1; i<val2 || (!(i<val2)&&!(val2<i)); i++) { // Use only the < operator of the template-type T.
                    srange.insert(i);
                }
            } else {
                return false;
            }
        } while (!str.empty());
    } catch (boost::bad_lexical_cast &e) {
        return false;
    }
    range.reserve(srange.size());
    for (std::set<std::size_t>::const_iterator it=srange.begin(); it!=srange.end(); it++) {
        range.push_back(*it);
    }
    return true;
}




/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename T, typename TIDX>
void tom::tools::subset(const std::vector<T> &v, const std::vector<TIDX> &vidx, std::vector<T> &vsub, bool unique, bool sort) {

    std::vector<TIDX> vsorted;
    const std::vector<TIDX> *pvidx = &vidx;
    if (unique) {
        std::set<TIDX> vidx_sorted;
        if (sort) {
            vidx_sorted.insert(vidx.begin(), vidx.end());
            vsorted.assign(vidx_sorted.begin(), vidx_sorted.end());
        } else {
            vsorted.reserve(vidx.size());
            for (typename std::vector<TIDX>::const_iterator it=vidx.begin(); it!=vidx.end(); it++) {
                if (vidx_sorted.insert(*it).second) {
                    vsorted.push_back(*it);
                }
            }
        }
        pvidx = &vsorted;
    } else if (sort) {
        std::multiset<TIDX> vidx_sorted(vidx.begin(), vidx.end());
        vsorted.assign(vidx_sorted.begin(), vidx_sorted.end());
        pvidx = &vsorted;
    }

    vsub.resize(0);
    vsub.reserve(pvidx->size());
    for (typename std::vector<TIDX>::const_iterator it=pvidx->begin(); it!=pvidx->end(); it++) {
        vsub.push_back(v[*it]);
    }
}


/****************************************************************************//**
 *
 *
 *******************************************************************************/
std::string tom::tools::gethostname() {
    #ifndef HOST_NAME_MAX
    #define HOST_NAME_MAX 256
    #endif
    std::vector<char> cbuff(HOST_NAME_MAX+2, 0);
    int i = ::gethostname(&cbuff[0], cbuff.size()-1);
    return i==0 ? &cbuff[0] : "";
}


/****************************************************************************//**
 *
 *
 *******************************************************************************/
template<typename T, typename _Compare>
void tom::tools::unify_shared_vector(std::vector<boost::shared_ptr<T> > &v, const _Compare &c) {

    std::map<T *, std::vector<std::size_t>, _Compare> m(c);

    typename std::vector<boost::shared_ptr<T> >::const_iterator vit;
    typename std::vector<std::size_t>::const_iterator vit2;
    typename std::map<T *, std::vector<std::size_t>, _Compare>::const_iterator mit;
    std::size_t i;

    for (i=0, vit=v.begin(); vit!=v.end(); vit++, i++) {
        std::vector<std::size_t> &vv = m[vit->get()];
        vv.push_back(i);
    }

    for (mit=m.begin(); mit!=m.end(); mit++) {
        const std::vector<std::size_t> &vv = mit->second;
        boost::shared_ptr<T> &p = v[vv[0]];
        for (vit2=vv.begin()+1; vit2!=vv.end(); vit2++) {
            v[*vit2] = p;
        }
    }
}




/***************************************************************************//**
 * Converts a floating point to a string with a certain format and precision.
 *
 * If the integer part is too large to fit into \c num_int, the precision of
 * the string will not be reduced, instead the string is longer than
 * num_int+1+num_prec characters.
 ******************************************************************************/
std::string tom::tools::align_decimal(double x, std::size_t num_int, std::size_t num_prec) {

    std::ostringstream ss;
    if (num_prec == 0) {
        ss << std::setw(num_int) << std::setprecision(0) << std::fixed << round(x);
        return ss.str();
    }


    ss << (boost::format("%_" + boost::lexical_cast<std::string>(num_int+1+num_prec) + "." + boost::lexical_cast<std::string>(num_prec) + "f") % x);

    std::string s = ss.str();

    // Replace trailing zeros...
    std::size_t i = s.find_last_not_of('0');
    if (i!=std::string::npos && i<s.size()) {
        if (s[i] == '.') {
            i--;
        }
        std::size_t l = s.size()-i-1;
        s.replace(i+1, l, l, ' ');
    }

    // Move the '-' before the number.
    if (s[0] == '-') {
        i = 1;
        while (i < s.size() && s[i]==' ') {
            i++;
        }
        if (i>1 && s[i-1]==' ') {
            s[0] = ' ';
            s[i-1] = '-';
        }
    }

    i = s.size();
    while (i > num_int+1+num_prec && s[i-1]==' ') {
        i--;
    }
    if (i < s.size()) {
        s.erase(i);
    }

    return s;

}




// template instantiations.
template bool tom::tools::str2range(const std::string &instr, std::vector<std::size_t> &range, char delimiter, char delimiter_range);


#include <tom/volume.hpp>
template void tom::tools::unify_shared_vector<tom::Volume<double>, tom::volume_less<double> >(std::vector<boost::shared_ptr<tom::Volume<double> > > &v, const tom::volume_less<double> &c);



template void tom::tools::subset<std::string, std::size_t>(const std::vector<std::string> &v, const std::vector<std::size_t> &vidx, std::vector<std::string> &vsub, bool unique, bool sort);

