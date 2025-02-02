#!/bin/sh

# Check if DISPLAY or WAYLAND_DISPLAY environment variables are set
if [ -z "$DISPLAY" ] && [ -z "$WAYLAND_DISPLAY" ]; then
  echo "No graphical display server is running."
  exit 1
fi

if [ -z "$(which evince)" ]; then
  echo "'evince' is not installed."
  exit 1
fi

# Check if 'evince' is installed
evince_installed=$(which evince)

if [ -z "$evince_installed" ]; then
  echo "'evince' is not installed."
  exit 1
fi

# Check if any of 'zathura', 'evince', or 'xpdf' is running with the same file
running_viewer=$(ps fx | grep -E 'zathura|evince|xpdf' | grep "$1" | grep -v grep)

# If no viewer is running with the file, start 'evince <filename>'
if [ -z "$running_viewer" ]; then
  command="evince \"$1\" 2>/dev/null &"
  echo "Running viewer: $command"
  eval $command
else
  echo "Viewer already running with file. Bringing it to focus."

  # Extract the PID of the running viewer
  viewer_pid=$(echo "$running_viewer" | awk '{print $1}')

  # Try to bring the viewer window to focus using xdotool, if installed
  if command -v xdotool > /dev/null; then
    xdotool windowactivate "$(xdotool search --pid $viewer_pid | tail -n 1)"
  else
    echo "xdotool is not installed; unable to focus window."
  fi
fi

