function expandAll() { 
  var wins = document.getElementsByTagName("DIV");
  for (var i = 0; i < wins.length; i++) { 
    if (wins[i].className == "cloak") { 
      wins[i].style.visibility = "visible";
      wins[i].style.display = "block";
    }
  }
} 

function collapseAll() { 
  var wins = document.getElementsByTagName("DIV");
  for (var i = 0; i < wins.length; i++) { 
    if (wins[i].className == "cloak") { 
      wins[i].style.visibility = "hidden";
      wins[i].style.display = "none";
    }
  }
} 

function showDiv(what) { 
  var wins = document.getElementsByTagName("DIV");
  for (var i = 0; i < wins.length; i++) { 
    if (wins[i].className == "cloak") { 
      if (wins[i].id == what) { 
        wins[i].style.visibility = "visible";
        wins[i].style.display = "block";
        current_menu = wins[i].id;
      }
    }
  }
} 

function toggleDiv(what) { 
  var wins = document.getElementsByTagName("DIV");
  for (var i = 0; i < wins.length; i++) { 
    if (wins[i].className == "cloak") { 
      if (wins[i].id == what) { 
        if (wins[i].style.visibility == "visible") { 
          wins[i].style.visibility = "hidden";
          wins[i].style.display = "none";
	} else { 
          wins[i].style.visibility = "visible";
          wins[i].style.display = "block";
	}
        current_menu = wins[i].id;
      }
    }
  }
} 
