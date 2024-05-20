// Initialize medium zoom.
$(document).ready(function () {
  var margin = window.innerWidth >= 992 ? window.innerWidth * 0.05 : 0;

  medium_zoom = mediumZoom("[data-zoomable]", {
    background: getComputedStyle(document.documentElement).getPropertyValue("--global-bg-color") + "ee", // + 'ee' for transparency.
    margin: margin,
  });

  // Update margin when window is resized.
  $(window).resize(function () {
    var margin = window.innerWidth >= 992 ? window.innerWidth * 0.05 : 0;
    medium_zoom.update({
      margin: margin,
    });
  });
});
