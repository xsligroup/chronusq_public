/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */

// Basically this timer is a graph of nodes that have up children separated
// by thread. All children only have one parent.
//
// A known issue with this design is that it does not allow for two overlapping
// timed regions with the same nested region - you will need to assign it to 
// one or the other.
//
// For example, if the execution graph looked like this:
//
//    Call level:  0   1
//                 --------- foo tick
//                 |
//                 |   ----- subfoo tick
//                 |   |
//                 |   ----- subfoo tock
//                 |
//                 |-------- bar tick
//                 |
//                 |   ----- baz tick
//                 |   |
//                 |   ----- baz tock
//                 |
//                 |-------- foo tock
//                 |
//                 --------- bar tock
//
// You would need to assign the baz subblock to the foo or the bar context,
// even though it is really part of both.
//
// In quantum chemistry programs, there aren't many reasons to have overlapping
// timing regions, so this issue shouldn't come up often.
//
// This timer also doesn't account for overhead. Experimentally, the overhead is
// in the tens of microseconds per tick/tock pair, but the exact number depends
// on your system and software environmental factors. In short, if you are
// optimizing/measuring something that takes microseconds and the overhead is
// substantial compared to your run time, you should use a timer more optimized
// for that case.
//
// Right now this timer can't handle nested OpenMP parallelism or other
// versions of SMP parallelism, but since the thread ids don't need to be
// sequential in each section, it can be extended to handle those cases.

#pragma once

#include <algorithm>
#include <cassert>
#include <chrono>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include <omp.h>

#include <util/singleton.hpp>


namespace ChronusQ {

  // Old time.hpp definitions - may want to phase out
#ifndef CQ_ENABLE_MPI
  using time_point = std::chrono::high_resolution_clock::time_point;
#else
  using time_point = double;
#endif


  static inline time_point tick() {

#ifndef CQ_ENABLE_MPI
    return std::chrono::high_resolution_clock::now();
#else
    return MPI_Wtime();
#endif

  }

  static inline double tock(const time_point& pt) {

    time_point now = tick();

#ifndef CQ_ENABLE_MPI
    return std::chrono::duration<double>(now - pt).count();
#else
    return (now - pt);
#endif

  }



  using CQMillisecond = std::chrono::duration<double, std::milli>;
  using CQSecond = std::chrono::duration<double>;
  using CQMinute = std::chrono::duration<double, std::ratio<60,1>>;

  // Forward declare
  class Timer;

  enum TimerUnit {
    MILLISECONDS,
    SECONDS,
    MINUTES
  };

  struct TimerOpts {
    TimerUnit unit{SECONDS};
    size_t dbgPrint{0};
    bool doSummary{true};
  };
  
  // Threadsafe wrapper around an unordered_map that also offers threadsafe
  // Id incrementing. 
  //
  // 0 is treated as a sentinel value by the Timer class, so incrementing
  //   starts at 1.
  template <typename IdT, typename ValT>
  class Registry {
  
    IdT nextId_ = 1;
    std::unordered_map<IdT,ValT> registry_;

    std::mutex id_mutex_;
    std::mutex reg_mutex_;
  
    public:
  
    IdT getId() {
      const std::lock_guard<std::mutex> lock(id_mutex_);
      IdT temp = nextId_++;
      return temp;
    }
  
    void add(IdT id, ValT val) {
      const std::lock_guard<std::mutex> lock(reg_mutex_);
      registry_.insert({id, val});
    }
  
    ValT operator[](IdT id) {
      const std::lock_guard<std::mutex> lock(reg_mutex_);
      ValT retval = registry_[id];
      return retval;
    };
  
  };
  
  // Basic node that represents a slice of time with a set of threaded sections
  // underneath
  class TimeSection {
  
    using IdType    = size_t;
    using Clock     = std::chrono::steady_clock;
    using TimePoint = std::chrono::time_point<Clock>;
    using RegType   = Registry<IdType,TimeSection*>;
  
    // Human readable label for this section
    std::string label_;
    // Unique Id for this section
    IdType id_;
  
    // Threads (outside) and subsections (inside)
    // The value is a unique_ptr so references to the object don't get
    //   invalidated on reallocation of the unordered map. 
    // Multiple shared_ptr would cause cyclical references, so we're not going to
    // do that
    std::unordered_map<size_t,std::unordered_map<IdType,std::unique_ptr<TimeSection>>> children_;
    TimeSection* parent_;
  
    // Registry to get ids and add subsections
    RegType* registry_;
  
    // Wall times for this section
    TimePoint start_;
    TimePoint stop_;

    bool stopped_ = false;

    std::mutex mutex_;
  
    friend class Timer;
  
    public:

    TimeSection(): label_(""), id_(0), parent_(nullptr), registry_(nullptr) {
      start_ = Clock::now();
    };
  
    TimeSection(
        std::string label,
        IdType id,
        RegType* reg,
        TimeSection* parent) :
      label_(label),
      id_(id),
      registry_(reg),
      parent_(parent)
    { 
      start_ = Clock::now();
    }

    //
    // Timing control
    //
  
    // Stop this section
    void stop() {
  
      if (stopped_) {
        return;
      }
  
      // Stop all subsections
      // Locked because tick can cause reallocation of children_
      //   on another thread - this is _very_ unlikely to happen during
      //   normal usage, but is possible.
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        for ( auto& thread: children_ )
          for ( auto& section: thread.second )
            section.second->stop();
      }
  
      stop_ = Clock::now();
      stopped_ = true;
  
    }
  
    // Start a new subsection
    IdType tick(std::string label, size_t threadId = 0) {
  
      IdType id = registry_->getId();
  
      auto child = std::make_unique<TimeSection>(label, id, registry_, this);
  
      registry_->add(id, child.get());

      // Locked due to insertion
      {
        const std::lock_guard<std::mutex> lock(mutex_);
        if (children_.find(threadId) == children_.end()) {
          std::unordered_map<IdType,std::unique_ptr<TimeSection>> temp;
          children_.insert(std::move(std::make_pair(threadId, std::move(temp))));
        }
        children_[threadId].insert(std::move(std::make_pair(id, std::move(child))));
      }
  
      return id;
    }
  
    // Stop a subsection
    void tock(size_t id, size_t threadId = 0) {
      // Locked because tick can cause reallocation of children_ on another
      // thread
      const std::lock_guard<std::mutex> lock(mutex_);
      children_[threadId].find(id)->second->stop();
    }
  
    //
    // Reporting
    //
    
    template <typename DurT>
    DurT duration() {
      return std::chrono::duration_cast<DurT>(stop_ - start_);
    }
  
    //
    // Searching ancestors (no recursion)
    //
  
    bool isAncestor(IdType id) {
      TimeSection* curr = this;
      while ( curr->parent_ != nullptr ) {
        curr = curr->parent_;
        if ( curr->id_ == id )
          return true;
      }
      return false;
    };
  
    TimeSection* findAncestor(std::string label) {
      TimeSection* curr = this;
      while ( curr->parent_ != nullptr ) {
        curr = curr->parent_;
        if ( curr->label_ == label )
          return curr;
      }
      return nullptr;
    };
  
  };
  
  //
  // Main Timer class
  //
  class Timer {
    public:

    using IdType = TimeSection::IdType;

    private:
  
    // Thread independent contexts
    std::vector<TimeSection*> contexts_;
    // Thread independent call levels
    std::vector<size_t> levels_;
  
    // Root of the tree
    TimeSection base_;
  
    // Registry of id to references
    Registry<IdType,TimeSection*> registry_;
  
    public:

    TimerOpts options;
  
    Timer() = default; 

    // Must call before use. Can't be in constructor because Timer must be
    //   default constructable (for Singleton ProgramTimer)
    void initialize(std::string baseLabel, size_t nThreads) {
      base_.label_ = baseLabel;
      base_.id_ = registry_.getId();
      base_.registry_ = &registry_;

      registry_.add(base_.id_, &base_);

      contexts_.resize(nThreads, nullptr);
      levels_.resize(nThreads, 1);

      contexts_[0] = &base_;
    };
  
    Timer(const Timer&) = delete;
  
    //
    // Section manipulation
    //
  
    // Start a new section
    IdType tick(std::string label, IdType parentId = 0, bool addContext = true) {
  
      size_t threadId = omp_get_thread_num();
  
      TimeSection* section = contexts_[threadId];
      if (parentId != 0)
        section = registry_[parentId];

      // Increment call level
      levels_[threadId]++;
  
      // Add the new section
      auto newId = section->tick(label, threadId);
  
      if (addContext) 
        contexts_[threadId] = registry_[newId];

      return newId;
  
    };
  
    // Stop a section by id
    void tock(IdType id = 0) {
  
      size_t threadId = omp_get_thread_num();
  
      TimeSection* section = contexts_[threadId];
      if (id != 0)
        section = registry_[id];
      section->stop();
  
      size_t s_id = section->id_;
  
      if ( section->parent_ != nullptr )
        if (contexts_[threadId]->id_ == s_id or
            contexts_[threadId]->isAncestor(s_id))
          contexts_[threadId] = registry_[section->parent_->id_];

      // Possibly print
      // This is bad for performance, but it is useful for debugging
      if ( levels_[threadId] < options.dbgPrint ) {
        std::stringstream message;
        std::string name = section->label_ + " duration: ";
        message << "  || " << std::left << std::setw(40) << name;
        message << std::right << std::setw(14) << std::fixed << std::setprecision(4);
        message << section->duration<CQMillisecond>().count() << " ms ||\n";
        std::cout << message.str();
      }

      // Decrement call level
      levels_[threadId]--;
  
    };
  
    // Stop a section by the best guess at which label it's referring to
    //
    // If current context matches label, stop that section, otherwise find the
    //   first section that matches the label in the history
    void tock(std::string label) {
  
      size_t threadId = omp_get_thread_num();
  
      TimeSection* section = contexts_[threadId];
      if ( section->label_ != label ) {
        section = section->findAncestor(label);
        assert(section != nullptr);
      }
  
      tock(section->id_);
    };


    // Time a functor with the given label
    template <typename F, typename... Args>
    auto timeOp(std::string label, const F& func, Args&&... args) {
      tick(label);
      func(std::forward<Args>(args)...);
      tock(label);
    }
  
    // 
    // Direct manipulation of current context
    //
    // XXX: These methods can make your history confusing/messed up if you aren't
    //   careful. If you choose to use these, you should keep track of the Ids
    //   returned by tick.
    //
    
    // Get FIRST instance of label's Id, starting from base
    IdType getLabelId(std::string label, IdType base = 0) {

      std::queue<TimeSection*> working;

      std::function<IdType(std::queue<TimeSection*>&)> tryNextId;
      tryNextId = [&](std::queue<TimeSection*>& next) -> IdType {
        TimeSection* curr = next.front();
        next.pop();
        if ( curr->label_ == label )
          return curr->id_;
        
        for ( auto& thread: curr->children_ )
        for ( auto& child: thread.second ) {
          next.push(child.second.get());
        }

        return 0;
      };

      size_t threadId = omp_get_thread_num();
      TimeSection* section = contexts_[threadId];
      if (base != 0)
        section = registry_[base];

      working.push(section);
      IdType id = 0;

      while ( !working.empty() ) {
        id = tryNextId(working);
        if ( id != 0 )
          break;
      }

      return id;

    };
    
    // Get Id of thread's current context
    IdType getContextId() {
      size_t threadId = omp_get_thread_num();
  
      assert(contexts_[threadId] != nullptr);
      return contexts_[threadId]->id_;
    };

    size_t getCallLevel() {
      return levels_[omp_get_thread_num()];
    };
  
    // Set context arbitrarily
    inline void setContext(IdType id, size_t level) {
      size_t threadId = omp_get_thread_num();
      contexts_[threadId] = registry_[id];
      levels_[threadId] = level;
    };
  
    // Return context to parent
    IdType popContext() {
      size_t threadId = omp_get_thread_num();
      IdType id = contexts_[threadId]->id_;
      contexts_[threadId] = contexts_[threadId]->parent_;
      return id;
    };
  
    //
    // Reporting
    //
  
    // Get duration of single section
    template <typename DurT>
    DurT getIdDuration(IdType id) {
      return registry_[id]->duration<DurT>();
    };
  
    // Get durations of all sections that match label
    template <typename DurT>
    std::vector<DurT> getDurations(std::string label, IdType base = 0,
      std::string filter = "") {

      auto tDurs = getThreadDurations<DurT>(label, base, filter);
  
      std::vector<DurT> results;
      for (auto& thread: tDurs)
        for (auto& dur: thread)
          results.push_back(dur);
  
      return results;
  
    };
  
    // Get (total, average) of all sections that match label
    template <typename DurT>
    std::pair<DurT,DurT> getDurationSummary(std::string label, IdType base = 0,
      std::string filter = "") {

      auto durations = getDurations<DurT>(label, base, filter);

      if (durations.size() == 0)
        return {DurT(0.), DurT(0.)};

      DurT sum(0);
      for (auto& dur: durations)
        sum += dur;

      return {sum, sum/durations.size()};
    }

    template <typename DurT>
    std::pair<DurT,DurT> getAggDurSummary(std::string label,
      std::string parentLabel, IdType base = 0) {

      auto durations = getAggregateDurations<DurT>(label, parentLabel, base);

      if (durations.size() == 0)
        return {DurT(0.), DurT(0.)};

      DurT sum(0.);
      for (auto& section: durations) {
        sum += *std::max_element(section.begin(), section.end());
      }

      return {sum, sum/durations.size()};

    }
  
    // Get average of all sections that match label
    template <typename DurT>
    DurT getDurationAvg(std::string label, IdType base = 0,
      std::string filter = "") {

      auto summary = getDurationSummary<DurT>(label, base, filter);
      return summary.second;
    }
  
    // Get total of all sections that match label
    template <typename DurT>
    DurT getDurationTotal(std::string label, IdType base = 0,
      std::string filter = "") {

      auto summary = getDurationSummary<DurT>(label, base, filter);
      return summary.first;
    }
  
    // Main get duration method
    //
    // If filterLabel is not "", it will only add durations that have an
    //   ancestor with that label
    template <typename DurT>
    std::vector<std::vector<DurT>> getThreadDurations(std::string label,
      IdType base = 0, std::string filterLabel = "") {

      std::function<void(size_t, size_t, bool, TimeSection*,
                         std::vector<DurT>&)> getTDurRecur;
      getTDurRecur =
        [&](size_t select, size_t tId, bool addDur, TimeSection* sect,
            std::vector<DurT>& results) {
    
          if (sect->label_ == label && select == tId && addDur) {
            results.push_back(sect->duration<DurT>());
          }

          addDur |= sect->label_ == filterLabel;
    
          for ( auto& thread: sect->children_ ) {
    
            if ( thread.first != select )
              continue;
    
            for ( auto& child: thread.second ) {
              getTDurRecur(select, thread.first, addDur, child.second.get(),
                           results);
            }
    
          }
    
        };
    
      std::vector<std::vector<DurT>> results;
    
      size_t threadId = omp_get_thread_num();
      TimeSection* sect = contexts_[threadId];
      if ( base != 0 )
        sect = registry_[base];
    
      for ( auto iTh = 0; iTh < contexts_.size(); iTh++ ) {
        std::vector<DurT> tempRes;

        // Handle the fact that we don't know which thread this section came 
        //   from by putting it on the calling thread stack
        getTDurRecur(iTh, threadId, filterLabel == "", sect, tempRes);

        results.push_back(std::move(tempRes));
      }
    
      return results;
    
    }

    // Aggregate the durations directly underneath parentLabel that match
    //   label.
    //
    // This can be useful in highly looped regions where you want to know
    //   how much time is spent in different parts of the loops but want to
    //   separate it by entrances to the loop (e.g. direct J/K builds:
    //   integral formation vs. density contraction)
    //
    // This does NOT give threads in the right order, so result[section[0]] may
    //   not correspond to thread 0.
    template <typename DurT>
    std::vector<std::vector<DurT>> getAggregateDurations(std::string label,
      std::string parentLabel, IdType base = 0) {

      std::function<void(TimeSection*,
                         std::vector<std::vector<DurT>>&)> getAggDurRecur;
      getAggDurRecur =
        [&](TimeSection* sect, std::vector<std::vector<DurT>>& results) {

          bool doAgg = sect->label_ == parentLabel;
          std::vector<DurT> current;

          for ( auto& thread: sect->children_ ) {
            DurT aggDur(0.);
            size_t tId = thread.first;
            for ( auto& child: thread.second ) {
              if ( doAgg && child.second->label_ == label ) {
                aggDur += child.second->duration<DurT>();
              }
              getAggDurRecur(child.second.get(), results);
            }
            current.push_back(aggDur);
          }

          if ( doAgg )
            results.push_back(current);

        };

      std::vector<std::vector<DurT>> results;

      size_t threadId = omp_get_thread_num();
      TimeSection* sect = contexts_[threadId];
      if ( base != 0 )
        sect = registry_[base];

      getAggDurRecur(sect, results);
      return results;

    }

    //
    // Printing
    //
  
    template <typename DurT>
    void printTimings(std::ostream& out, std::string unit = "", IdType id = 0,
      int maxDepth = 3) {
  
      TimeSection* printBase = contexts_[0];
      if (id != 0)
        printBase = registry_[id];
      auto total = printBase->duration<DurT>();
  
      // Print header
      std::string first = "Section Label";
      std::string second = "  Duration (";
      second += unit;
      second += ")  Percent";
  
      // Use a "max rational section label length" of 20. If it's longer than
      // that, the printing won't be pretty but will still work
      size_t fieldWidth = 20;
      size_t width = maxDepth*2 + fieldWidth;
  
      // Size of the data fields
      size_t fixed = 22;
      size_t dataw = std::max(fixed, second.size());
  
      // Size of max width
      size_t outw = width + dataw;
  
      //
      // Recursive printing definition
      //
      std::function<void(TimeSection*, size_t)> recursivePrint;
      recursivePrint = [&](TimeSection* sect, size_t level) {
        if ( level > maxDepth )
          return;
  
        DurT dur = sect->duration<DurT>();
  
        // Print this section
        std::stringstream ss;
        for ( auto i = 0; i < level; i++ )
          ss << "| ";
        ss << sect->label_ << " ";
        out << std::setw(width);
        out << std::left << std::setfill('-') << ss.str();
  
        out << std::right << std::setfill(' ');
        out << std::setw(outw - width - fixed) << "";
  
        out << ' ' << std::setw(12) << std::fixed << std::setprecision(4);
        out << dur.count();
  
        out << "  " << std::setw(7) << std::fixed << std::setprecision(3);
        out << dur.count()/total.count()*100;
  
        out << std::endl;
  
  
        auto sortByStart = [](TimeSection* a, TimeSection* b) {
          return (a->start_ < b->start_);
        };
  
        // Recurse into children
        bool threaded = sect->children_.size() > 1;
  
        std::vector<TimeSection*> children;
  
        for ( auto iTh = 0; iTh < contexts_.size(); iTh++ ) {
  
          if (sect->children_.find(iTh) == sect->children_.end())
            continue;
  
          children.clear();
  
          for (auto& child: sect->children_[iTh] ) {
            children.push_back(child.second.get());
          }
  
          sort(children.begin(), children.end(), sortByStart);
  
          if (threaded) {
            std::stringstream ss;
            for ( auto i = 0; i < level+1; i++ )
              ss << "| ";
            ss << "*** Thread " << iTh << " ";
            out << std::left << std::setfill('*');
            out << std::setw(outw) << ss.str() << std::endl;
          }
  
          for ( auto child: children ) {
            recursivePrint(child, level+1);
          }
        }
      }; // Recursive printing definition
  
      out << std::setw(width) << std::left << std::setfill(' ') << first;
      out << std::right << std::setw(dataw-10) << second << std::endl;
      out << std::setw(outw) << std::setfill('=') << "" << std::endl;
      recursivePrint(printBase, 0);
      out << std::setw(outw) << std::setfill('=') << "" << std::endl;
      out << std::setfill(' ');
  
    } // printTimings

    void printTimings(std::ostream& out, IdType id = 0) {
      if ( options.unit == MILLISECONDS )
        printTimings<CQMillisecond>(out, "ms", id, options.dbgPrint);
      else if ( options.unit == SECONDS )
        printTimings<CQSecond>(out, "s", id, options.dbgPrint);
      else if ( options.unit == MINUTES )
        printTimings<CQMinute>(out, "min", id, options.dbgPrint);
    };
  
  }; // class Timer

  // Singleton for the program-wide timer
  class ProgramTimer: public Singleton<Timer> {

    using IdType = Timer::IdType;

    public:

    static void initialize(std::string tag, size_t nThreads) {
      instance()->initialize(tag, nThreads);
    }

    static IdType tick(std::string tag, IdType parentId = 0,
      bool addContext = true) {
      return instance()->tick(tag, parentId, addContext);
    }

    static void tock(IdType id = 0) {
      instance()->tock(id);
    }

    static void tock(std::string tag) {
      instance()->tock(tag);
    }

    template <typename F, typename... Args>
    static auto timeOp(std::string tag, const F& func, Args&&... args) {
      return instance()->timeOp(tag, func, std::forward<Args>(args)...);
    }

    static void setContext(IdType id, size_t level) {
      instance()->setContext(id, level);
    }

    static IdType getLabelId(std::string label, IdType base = 0) {
      return instance()->getLabelId(label, base);
    }

    static IdType getContextId() {
      return instance()->getContextId();
    }

    static size_t getCallLevel() {
      return instance()->getCallLevel();
    }

    template <typename DurT>
    static DurT getIdDuration(IdType id) {
      return instance()->getIdDuration<DurT>(id);
    }

    template <typename DurT>
    static std::vector<DurT> getDurations(std::string tag, IdType base = 0,
      std::string filter = "") {
      return instance()->getDurations<DurT>(tag, base, filter);
    }

    template <typename DurT>
    static std::pair<DurT,DurT> getDurationSummary(std::string tag,
      IdType base = 0, std::string filter = "") {
      return instance()->getDurationSummary<DurT>(tag, base, filter);
    }

    template <typename DurT>
    static DurT getDurationTotal(std::string tag, IdType base = 0,
      std::string filter = "") {
      return instance()->getDurationTotal<DurT>(tag, base, filter);
    }

    template <typename DurT>
    static DurT getDurationAvg(std::string tag, IdType base = 0,
      std::string filter = "") {
      return instance()->getDurationAvg<DurT>(tag, base, filter);
    }

    template <typename DurT>
    static std::vector<std::vector<DurT>> getThreadDurations(std::string tag,
      IdType base = 0, std::string filter = "") {
      return instance()->getDurationAvg<DurT>(tag, base, filter);
    }

    template <typename DurT>
    static std::vector<std::vector<DurT>> getAggregateDurations(std::string tag,
      std::string parentTag, IdType base = 0) {
      return instance()->getDurationAvg<DurT>(tag, parentTag, base);
    }

    template <typename DurT>
    static std::pair<DurT,DurT> getAggDurSummary(std::string tag,
      std::string parentTag, IdType base = 0) {
      return instance()->getAggDurSummary<DurT>(tag, parentTag, base);
    }

    template <typename DurT>
    static void printTimings(std::ostream& out, std::string unit = "",
      IdType id = 0, int maxDepth = 3) {
      instance()->printTimings<DurT>(out, unit, id, maxDepth);
    }

    static void printTimings(std::ostream& out, IdType id = 0) {
      instance()->printTimings(out, id);
    }

  };

} // namespace ChronusQ
