// Copyright (C) 2005 Universitat d'Alacant / Universidad de Alicante
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, see <http://www.gnu.org/licenses/>.

#ifndef EXCEPTION_APERTIUM_TAGGER_H
#define EXCEPTION_APERTIUM_TAGGER_H

#include "exception_type.h"

#include <sstream>

namespace Apertium {
namespace Exception {

#define EXCEPTION(EXCEPTION_TYPE)                                              \
  class EXCEPTION_TYPE : public ::Apertium::ExceptionType {                    \
  public:                                                                      \
    EXCEPTION_TYPE(const char *const what_) : ExceptionType(what_) {}          \
    EXCEPTION_TYPE(const std::string &what_) : ExceptionType(what_) {}         \
    EXCEPTION_TYPE(const std::stringstream &what_) : ExceptionType(what_) {}   \
    ~EXCEPTION_TYPE() throw() {}                                               \
  };

namespace atools {
EXCEPTION(UnexpectedOption)
}

namespace Optional {
EXCEPTION(TheOptionalTypePointer_null)
}

#undef EXCEPTION
}
}

#endif // EXCEPTION_APERTIUM_TAGGER_H
