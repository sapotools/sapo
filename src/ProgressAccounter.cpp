/**
 * @file ProgressAccounter.cpp
 * This file declares a class to account computation progresses.
 *
 * @author Alberto Casagrande <acasagrande@units.it>
 * @version 0.1
 */

#include "ProgressAccounter.h"

const unsigned int &
ProgressBar::increase_performed(const unsigned int delta_steps)
{
  std::unique_lock<std::mutex> lock(_mutex);

  return unsafe_increase_performed_to(_performed + delta_steps);
}

const unsigned int &
ProgressBar::increase_performed_to(const unsigned int performed)
{
  std::unique_lock<std::mutex> lock(_mutex);

  return unsafe_increase_performed_to(performed);
}

const unsigned int &
ProgressBar::unsafe_increase_performed_to(const unsigned int performed)
{
  if (performed < _performed || _represented_steps == _expected) {
    return _performed;
  }

  if (_represented_steps == 0 && _represented_steps != _expected) {
    _bstream << _preamble << "0% " << std::flush;
  }

  _performed = performed;

  if (_num_of_bar_dots > 0 && _expected > 0) {
    float dot_per_step = float(_num_of_bar_dots) / _expected;

    const unsigned int already_displayed
        = int(dot_per_step * _represented_steps);
    const unsigned int to_display = int(dot_per_step * _performed);

    for (unsigned int dot = already_displayed + 1; dot <= to_display; ++dot) {
      _bstream << "#" << std::flush;
      if (_num_of_bar_dots > 10
          && int((dot - 1) * 5.0 / _num_of_bar_dots)
                 < int(dot * 5.0 / _num_of_bar_dots)) {
        _bstream << " " << int(dot * 100.0 / _num_of_bar_dots) << "% "
                 << std::flush;
      }
    }
  }
  _represented_steps = _performed;

  if (_performed == _expected) {
    _bstream << std::endl << std::flush;
  }

  return _performed;
}