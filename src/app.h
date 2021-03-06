// taken from http://code.google.com/p/cpp-project-template/
#ifndef HEADER_SRC_APP_H_INCLUDED
#define HEADER_SRC_APP_H_INCLUDED

#include <QtCore>
#include <QtGui>

class App : public QApplication
{
        Q_OBJECT
    public:
        App(int &argc, char **argv);
        ~App();
        
        App *INSTANCE();

        QString getProjectName();
        QString getProjectCodeName();
        QString getProjectVendorID();
        QString getProjectVendorName();
        QString getProjectID();
        int getProjectMajorVersion();
        int getProjectMinorVersion();
        int getProjectPatchVersion();
        QString getProjectVersion();
        QString getProjectCopyrightYears();
        QString getProjectInvocation();

    private:
        void interactiveMain(const std::string& o_f_n);
        void consoleMain(const std::string& o_f_n);
    
        void printHelpMessage();
        void printVersionMessage();
        void printVersionTripletMessage();
        void printApplicationIdentifier();
        void setPreference(const std::string &key, const std::string &val);
        void unsetPreference(const std::string &key);
        void printPreference(const std::string &key)const;
        void printAllPreferences()const;
        std::string getKeyName(const std::string &key)const;
        std::string getKeyRepr(const std::string &key)const;
        std::string convert(const QString &str)const;
        QString convert(const std::string &str)const;
        
        
        static App *_instance;
        QString _invocation;
        bool _interactive;
};

#endif
