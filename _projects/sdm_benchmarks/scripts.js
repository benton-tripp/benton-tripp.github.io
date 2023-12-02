$(document).ready(function(){
  $('.shared-menu-item').click(function(e) {
    $('.shared-menu-item').removeClass('active');
    $(this).addClass('active');
    
    // Navigate to the target tab
    var tar_a = $(this).children()[0];
    var target = $(tar_a).attr('href');
    $('.tab-pane').removeClass('active');
    $(target).addClass('active');
  });
});