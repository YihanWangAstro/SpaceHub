#ifndef ERRHAND_H
#define ERRHAND_H

#include <iostream>
#include <cstring>
#include <sstream>
#include "macros.h"
#ifdef TELEGRAM
#include <bot_alert.h>
#endif
namespace NOTICE
{
#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define NEWLINE printf("\n");
    constexpr size_t WIDTH = 80;
    extern bool Message;
    inline void Telegram(const char* host, const char* msg)
    {
        if(Message)
            printf(ANSI_COLOR_RED  "%15s "  ANSI_COLOR_RESET ": %s\n", host, msg);
        
    };
    
    inline void Title(const char* T)
    {
        if(Message)
        {
            size_t len = strlen(T);
            size_t rest = (WIDTH - len) >> 1;
            for(size_t i = 0 ; i < rest; i++)printf("=");
            printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET, T);
            for(size_t i = 0 ; i < rest; i++)printf("=");
            NEWLINE
        }
    };
    
    inline void SubTitle(const char* T)
    {
        if(Message)
        {
            size_t len = strlen(T);
            size_t rest = (WIDTH - len) >> 1;
            for(size_t i = 0 ; i < rest; i++)printf("-");
            printf(ANSI_COLOR_RED "%s" ANSI_COLOR_RESET,T);
            for(size_t i = 0 ; i < rest; i++)printf("-");
            NEWLINE
        }
    };
    
    inline void EraseLine()
    {
        for(size_t i = 0 ; i < WIDTH; i++)printf("\b");
    };
    
    inline void Line()
    {
        if(Message)
        {
            for(size_t i = 0 ; i < WIDTH; i++)printf("=");
            NEWLINE
        }
    };
    
    inline void SubLine()
    {
        if(Message)
        {
            for(size_t i = 0 ; i < WIDTH; i++)printf("-");
            NEWLINE
        }
    };
    
    inline void RunInfo(double timeLimit, double outputsize_terval, double tolerance)
    {
        if(Message)
        {
            if(outputsize_terval != 0)
                printf("Time Limit = " ANSI_COLOR_GREEN "%-8.3le" ANSI_COLOR_RESET
                       " year    Output size_terval = " ANSI_COLOR_GREEN "%-8.3le" ANSI_COLOR_RESET
                       " year    Error tolerance = " ANSI_COLOR_GREEN "%-8.3le" ANSI_COLOR_RESET
                       "\n", timeLimit/YEAR, outputsize_terval/YEAR, tolerance);
            else
                printf("Time Limit = " ANSI_COLOR_GREEN "%-8.3le" ANSI_COLOR_RESET
                       " year    No output File    Error tolerance = " ANSI_COLOR_GREEN
                       "%-8.3le" ANSI_COLOR_RESET "\n", timeLimit/YEAR, tolerance);
        }
    };
}

class errhand
{
  public:
    errhand(std::string err_msg_input, const char* file_input, size_t line_input)
            {
                err_msg = err_msg_input;
                file    = to_string_loc(file_input);
                line    = line_input;
         
                print_to_stdout();
                invoke_telegram_bot();
            }
    std::string get_msg()  const {return err_msg;}
    std::string get_file() const {return file;   }
    size_t    get_line() const {return line;   }
    std::string to_string_loc(const char* obj)
            {
                std::ostringstream out;
                out << obj;
                return out.str();
            }
    
       void invoke_telegram_bot()
            {
#ifdef TELEGRAM
     
                bot = new bot_alert;
                string message = "<b>An error occurred</b>\n<b>Message</b> = " + get_msg() + "\n" + "<b>File</b> = "
                               + get_file() + "\n" + "<b>Line</b> = " + bot->to_string_loc(get_line()) + "\n";
         
                bot->init();
                bot->send_message(message);
                bot->send_sticker("sad");
#endif
                return;
            }

       void print_to_stdout()
            {
                std::cout << " An error occurred: " << '\n';
                std::cout << "    Message >>> " << get_msg() << '\n';
                std::cout << "    File    >>> " << get_file() << '\n';
                std::cout << "    Line    >>> " << get_line() << std::endl;
    
        return;
            }
  private:
    std::string   err_msg;
    std::string   file;
    size_t      line;
#ifdef TELEGRAM
    bot_alert* bot;
#endif

};



#endif
