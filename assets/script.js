$(document).ready(function () {

    // Toggle menu on hamburger click
    $('#navbar-toggle').on('click', function (e) {
        e.stopPropagation(); // Prevent click from bubbling up to document
        $('#navbar-menu').toggleClass('show');
    });

    // Close menu when clicking outside of it
    $(document).on('click', function (e) {
        // If the click isn’t inside #navbar-menu or #navbar-toggle, remove 'show'
        if (!$(e.target).closest('#navbar-menu, #navbar-toggle').length) {
            $('#navbar-menu').removeClass('show');
        }
    });

    // Also close the menu when clicking any .nav-item
    $('.navbar-item').on('click', function () {
        $('#navbar-menu').removeClass('show');
    });
});
