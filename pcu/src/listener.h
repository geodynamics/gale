/*
**  Copyright 2008 Luke Hodkinson
**
**  This file is part of pcu.
**
**  pcu is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  Foobar is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License for more details.
**
**  You should have received a copy of the GNU General Public License
**  along with Foobar.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef pcu_listener_h
#define pcu_listener_h

struct pcu_listener_t {
      pcu_suiteentry_t* suitebegin;
      pcu_suiteentry_t* suiteend;
      pcu_testentry_t* testbegin;
      pcu_testentry_t* testend;
      pcu_checkentry_t* checkdone;
      void* data;
};

#endif
