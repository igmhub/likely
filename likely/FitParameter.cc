// Created 28-May-2012 by David Kirkby (University of California, Irvine) <dkirkby@uci.edu>

#include "likely/FitParameter.h"
#include "likely/RuntimeError.h"

#include "boost/format.hpp"
#include "boost/function.hpp"
#include "boost/bind.hpp"
#include "boost/regex.hpp"
#include "boost/spirit/include/qi.hpp"
#include "boost/spirit/include/phoenix_core.hpp"
#include "boost/spirit/include/phoenix_operator.hpp"
#include "boost/spirit/include/phoenix_stl.hpp"

#include <iterator>
#include <iostream>

namespace local = likely;

namespace qi = boost::spirit::qi;
namespace ascii = boost::spirit::ascii;
namespace phoenix = boost::phoenix;

std::string const &local::FitParameter::getValidNameCharacters() {
    // The hyphen must appear as the last char below for correct interpretation.
    static std::string charset(" a-zA-Z0-9()*/+-");
    return charset;
}

local::FitParameter::FitParameter(std::string const &name, double value, double error)
: _name(name)
{
    boost::cmatch what;
    boost::regex pattern(std::string("[") + getValidNameCharacters() + "]+");
    if(!boost::regex_match(name,pattern)) {
        throw RuntimeError("FitParameter: name uses invalid characters \"" + name + "\"");
    }
    setValue(value);
    setError(error);
}

void local::FitParameter::setError(double error) {
    if(error < 0) {
        throw RuntimeError("FitParameter::setError: error must be >= 0.");
    }
    _error = error;
}

void local::getFitParameterValues(FitParameters const &parameters, Parameters &values, bool onlyFloating) {
    values.resize(0);
    if(!onlyFloating) values.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) values.push_back(iter->getValue());
    }
}

void local::getFitParameterErrors(FitParameters const &parameters, Parameters &errors, bool onlyFloating) {
    errors.resize(0);
    if(!onlyFloating) errors.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) errors.push_back(iter->getError());
    }
}

void local::getFitParameterNames(FitParameters const &parameters, std::vector<std::string> &names,
bool onlyFloating) {
    names.resize(0);
    if(!onlyFloating) names.reserve(parameters.size());
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(!onlyFloating || iter->isFloating()) names.push_back(iter->getName());
    }
}

int local::countFloatingFitParameters(FitParameters const &parameters) {
    int count(0);
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(iter->isFloating()) count++;
    }
    return count;
}

void local::printFitParametersToStream(FitParameters const &parameters, std::ostream &out,
std::string const &formatSpec) {
    boost::format formatter(formatSpec), label("%20s = ");
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        out << (label % iter->getName()) << (formatter % iter->getValue());
        if(iter->isFloating()) {
            out << " +/- " << formatter % iter->getError();
        }
        out << std::endl;
    }
}

int local::findFitParameterByName(FitParameters const &parameters, std::string const &name) {
    for(FitParameters::const_iterator iter = parameters.begin(); iter != parameters.end(); ++iter) {
        if(name == iter->getName()) return std::distance(parameters.begin(),iter);
    }
    throw RuntimeError("findFitParameterByName: name not found: '" + name + "'");
}

namespace likely {
namespace fitpar {
    // Declare our script grammar.
    struct Grammar : qi::grammar<std::string::const_iterator,ascii::space_type> {

        Grammar(FitParameters &_params) : params(_params), base_type(script) {

            using qi::double_;
            using qi::lit;
            using qi::no_skip;
            using qi::char_;

            script = command >> *(';' >> command);

            command =
                ( "value" >> name >> '=' >> double_[boost::bind(&Grammar::setValue,this,::_1)] ) |
                ( "error" >> name >> '=' >> double_[boost::bind(&Grammar::setError,this,::_1)] ) |
                ( "fix" >> name[boost::bind(&Grammar::fix,this)] ) |
                ( "release" >> name[boost::bind(&Grammar::release,this)] );
            
            name = lit('[')[boost::bind(&Grammar::beginName,this)]
                >> no_skip[+char_(FitParameter::getValidNameCharacters())[
                    boost::bind(&Grammar::addToName,this,::_1)]]
                >> lit(']')[boost::bind(&Grammar::endName,this)];
            
        }
        
        qi::rule<std::string::const_iterator, ascii::space_type> script, command, name;        

        FitParameters &params;
        std::string theName;
        int index;
        
        void beginName() {
            theName.clear();
        }
        void addToName(char c) {
            theName += c;
        }
        void endName() {
            index = findFitParameterByName(params,theName);
        }
        void setValue(double value) {
            params[index].setValue(value);
        }
        void setError(double error) {
            params[index].setError(error);
        }
        void fix() {
            params[index].fix();
        }
        void release() {
            params[index].release();
        }
    };
} // grammar
} // likely

void local::modifyFitParameters(FitParameters &parameters, std::string const &script) {
    FitParameters modified(parameters);
    fitpar::Grammar grammar(modified);
    std::string::const_iterator iter = script.begin();
    bool ok = qi::phrase_parse(iter, script.end(), grammar, ascii::space);
    if(!ok || iter != script.end()) {
        throw RuntimeError("modifyFitParameters: syntax error in script.");
    }
    parameters = modified;
}
