// modules are defined as an array
// [ module function, map of requires ]
//
// map of requires is short require name -> numeric require
//
// anything defined in a previous bundle is accessed via the
// orig method which is the require for previous bundles
parcelRequire = (function (modules, cache, entry, globalName) {
  // Save the require from previous bundle to this closure if any
  var previousRequire = typeof parcelRequire === 'function' && parcelRequire;
  var nodeRequire = typeof require === 'function' && require;

  function newRequire(name, jumped) {
    if (!cache[name]) {
      if (!modules[name]) {
        // if we cannot find the module within our internal map or
        // cache jump to the current global require ie. the last bundle
        // that was added to the page.
        var currentRequire = typeof parcelRequire === 'function' && parcelRequire;
        if (!jumped && currentRequire) {
          return currentRequire(name, true);
        }

        // If there are other bundles on this page the require from the
        // previous one is saved to 'previousRequire'. Repeat this as
        // many times as there are bundles until the module is found or
        // we exhaust the require chain.
        if (previousRequire) {
          return previousRequire(name, true);
        }

        // Try the node require function if it exists.
        if (nodeRequire && typeof name === 'string') {
          return nodeRequire(name);
        }

        var err = new Error('Cannot find module \'' + name + '\'');
        err.code = 'MODULE_NOT_FOUND';
        throw err;
      }

      localRequire.resolve = resolve;
      localRequire.cache = {};

      var module = cache[name] = new newRequire.Module(name);

      modules[name][0].call(module.exports, localRequire, module, module.exports, this);
    }

    return cache[name].exports;

    function localRequire(x){
      return newRequire(localRequire.resolve(x));
    }

    function resolve(x){
      return modules[name][1][x] || x;
    }
  }

  function Module(moduleName) {
    this.id = moduleName;
    this.bundle = newRequire;
    this.exports = {};
  }

  newRequire.isParcelRequire = true;
  newRequire.Module = Module;
  newRequire.modules = modules;
  newRequire.cache = cache;
  newRequire.parent = previousRequire;
  newRequire.register = function (id, exports) {
    modules[id] = [function (require, module) {
      module.exports = exports;
    }, {}];
  };

  var error;
  for (var i = 0; i < entry.length; i++) {
    try {
      newRequire(entry[i]);
    } catch (e) {
      // Save first error but execute all entries
      if (!error) {
        error = e;
      }
    }
  }

  if (entry.length) {
    // Expose entry point to Node, AMD or browser globals
    // Based on https://github.com/ForbesLindesay/umd/blob/master/template.js
    var mainExports = newRequire(entry[entry.length - 1]);

    // CommonJS
    if (typeof exports === "object" && typeof module !== "undefined") {
      module.exports = mainExports;

    // RequireJS
    } else if (typeof define === "function" && define.amd) {
     define(function () {
       return mainExports;
     });

    // <script>
    } else if (globalName) {
      this[globalName] = mainExports;
    }
  }

  // Override the current require with this new one
  parcelRequire = newRequire;

  if (error) {
    // throw error from earlier, _after updating parcelRequire_
    throw error;
  }

  return newRequire;
})({"node_modules/parcel/src/builtins/bundle-url.js":[function(require,module,exports) {
var bundleURL = null;

function getBundleURLCached() {
  if (!bundleURL) {
    bundleURL = getBundleURL();
  }

  return bundleURL;
}

function getBundleURL() {
  // Attempt to find the URL of the current script and use that as the base URL
  try {
    throw new Error();
  } catch (err) {
    var matches = ('' + err.stack).match(/(https?|file|ftp|chrome-extension|moz-extension):\/\/[^)\n]+/g);

    if (matches) {
      return getBaseURL(matches[0]);
    }
  }

  return '/';
}

function getBaseURL(url) {
  return ('' + url).replace(/^((?:https?|file|ftp|chrome-extension|moz-extension):\/\/.+)\/[^/]+$/, '$1') + '/';
}

exports.getBundleURL = getBundleURLCached;
exports.getBaseURL = getBaseURL;
},{}],"node_modules/parcel/src/builtins/css-loader.js":[function(require,module,exports) {
var bundle = require('./bundle-url');

function updateLink(link) {
  var newLink = link.cloneNode();

  newLink.onload = function () {
    link.remove();
  };

  newLink.href = link.href.split('?')[0] + '?' + Date.now();
  link.parentNode.insertBefore(newLink, link.nextSibling);
}

var cssTimeout = null;

function reloadCSS() {
  if (cssTimeout) {
    return;
  }

  cssTimeout = setTimeout(function () {
    var links = document.querySelectorAll('link[rel="stylesheet"]');

    for (var i = 0; i < links.length; i++) {
      if (bundle.getBaseURL(links[i].href) === bundle.getBundleURL()) {
        updateLink(links[i]);
      }
    }

    cssTimeout = null;
  }, 50);
}

module.exports = reloadCSS;
},{"./bundle-url":"node_modules/parcel/src/builtins/bundle-url.js"}],"node_modules/@fortawesome/fontawesome-free/css/all.css":[function(require,module,exports) {
var reloadCSS = require('_css_loader');

module.hot.dispose(reloadCSS);
module.hot.accept(reloadCSS);
},{"./..\\webfonts\\fa-brands-400.eot":[["fa-brands-400.1bb139e6.eot","node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.eot"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.eot"],"./..\\webfonts\\fa-brands-400.woff2":[["fa-brands-400.1d34615d.woff2","node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.woff2"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.woff2"],"./..\\webfonts\\fa-brands-400.woff":[["fa-brands-400.eca31406.woff","node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.woff"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.woff"],"./..\\webfonts\\fa-brands-400.ttf":[["fa-brands-400.df86de32.ttf","node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.ttf"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.ttf"],"./..\\webfonts\\fa-brands-400.svg":[["fa-brands-400.f1eb0e8c.svg","node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.svg"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-brands-400.svg"],"./..\\webfonts\\fa-regular-400.eot":[["fa-regular-400.a2c1909d.eot","node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.eot"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.eot"],"./..\\webfonts\\fa-regular-400.woff2":[["fa-regular-400.5ca8c932.woff2","node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.woff2"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.woff2"],"./..\\webfonts\\fa-regular-400.woff":[["fa-regular-400.3c3cc54e.woff","node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.woff"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.woff"],"./..\\webfonts\\fa-regular-400.ttf":[["fa-regular-400.cde05ce7.ttf","node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.ttf"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.ttf"],"./..\\webfonts\\fa-regular-400.svg":[["fa-regular-400.6ef294e6.svg","node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.svg"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-regular-400.svg"],"./..\\webfonts\\fa-solid-900.eot":[["fa-solid-900.90890cef.eot","node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.eot"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.eot"],"./..\\webfonts\\fa-solid-900.woff2":[["fa-solid-900.da0e0451.woff2","node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.woff2"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.woff2"],"./..\\webfonts\\fa-solid-900.woff":[["fa-solid-900.935b31ea.woff","node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.woff"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.woff"],"./..\\webfonts\\fa-solid-900.ttf":[["fa-solid-900.f2409036.ttf","node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.ttf"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.ttf"],"./..\\webfonts\\fa-solid-900.svg":[["fa-solid-900.c87ba59a.svg","node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.svg"],"node_modules/@fortawesome/fontawesome-free/webfonts/fa-solid-900.svg"],"_css_loader":"node_modules/parcel/src/builtins/css-loader.js"}],"Constants.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.Software = {
  name: 'IoT Workbench',
  abbreviation: 'IW',
  version: '0.0.1'
};
exports.Font = {
  label: '9px Arial',
  normal: '10px Arial',
  axisName: '12px Arial',
  highlight: '14px Arial'
};
},{}],"User.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var User =
/** @class */
function () {
  function User(firstName, middleInitial, lastName) {
    this.firstName = firstName;
    this.middleInitial = middleInitial;
    this.lastName = lastName;

    if (middleInitial) {
      this.fullName = firstName + " " + middleInitial + " " + lastName;
    } else {
      this.fullName = firstName + " " + lastName;
    }
  }

  return User;
}();

exports.User = User;
},{}],"Workbench.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Main_1 = require("./Main");

var Workbench =
/** @class */
function () {
  function Workbench(canvasId) {
    this.gridSize = 20;
    this.canvas = document.getElementById(canvasId);
    this.canvas.addEventListener("mousedown", this.mouseDown.bind(this), false);
    this.canvas.addEventListener("mouseup", this.mouseUp.bind(this), false);
    this.canvas.addEventListener("mousemove", this.mouseMove.bind(this), false);
    this.canvas.addEventListener('contextmenu', this.openContextMenu.bind(this), false);
  }

  Workbench.prototype.draw = function () {
    var context = this.canvas.getContext('2d');
    context.clearRect(0, 0, this.canvas.width, this.canvas.height);
    this.drawGrid(context);
  };

  Workbench.prototype.drawGrid = function (context) {
    context.beginPath();
    context.strokeStyle = "LightSkyBlue";

    for (var i = 1; i <= this.canvas.height / this.gridSize; i++) {
      context.moveTo(0, i * this.gridSize);
      context.lineTo(this.canvas.width, i * this.gridSize);
    }

    for (var i = 1; i <= this.canvas.width / this.gridSize; i++) {
      context.moveTo(i * this.gridSize, 0);
      context.lineTo(i * this.gridSize, this.canvas.height);
    }

    context.stroke();
    context.closePath();
    context.restore();
  }; // detect if (x, y) is inside this workbench


  Workbench.prototype.contains = function (x, y) {
    return x > this.canvas.offsetLeft && x < this.canvas.offsetLeft + this.canvas.width && y > this.canvas.offsetTop && y < this.canvas.offsetTop + this.canvas.height;
  };

  Workbench.prototype.getX = function () {
    return 10;
  };

  Workbench.prototype.getY = function () {
    return 10;
  };

  Workbench.prototype.getWidth = function () {
    return this.canvas.width;
  };

  Workbench.prototype.getHeight = function () {
    return this.canvas.height;
  };

  Workbench.prototype.mouseDown = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
  };

  Workbench.prototype.mouseUp = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
    var context = this.canvas.getContext("2d");
  };

  Workbench.prototype.mouseMove = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
    var context = this.canvas.getContext("2d");
  };

  Workbench.prototype.openContextMenu = function (e) {
    e.preventDefault();
    var menu = document.getElementById(Main_1.contextMenus.workbench.id);
    menu.style.left = e.clientX + "px";
    menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
    menu.classList.add("show-menu");
  };

  return Workbench;
}();

exports.Workbench = Workbench;
},{"./Main":"Main.ts"}],"components/Board.ts":[function(require,module,exports) {
"use strict";
/*
 * A rectangular printed circuit board (PCB)
 *
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Board =
/** @class */
function () {
  function Board(canvasId) {
    this.handles = [];
    this.canvas = document.getElementById(canvasId);
  }

  Board.prototype.getUid = function () {
    return this.uid;
  };

  Board.prototype.getX = function () {
    return this.canvas.offsetLeft;
  };

  Board.prototype.setX = function (x) {
    this.canvas.style.left = x + "px";
  };

  Board.prototype.getY = function () {
    return this.canvas.offsetTop;
  };

  Board.prototype.setY = function (y) {
    this.canvas.style.top = y + "px";
  };

  Board.prototype.getWidth = function () {
    return this.canvas.width;
  };

  Board.prototype.getHeight = function () {
    return this.canvas.height;
  }; // detect if (x, y) is inside this board


  Board.prototype.contains = function (x, y) {
    return x > this.canvas.offsetLeft && x < this.canvas.offsetLeft + this.canvas.width && y > this.canvas.offsetTop && y < this.canvas.offsetTop + this.canvas.height;
  };

  Board.prototype.whichHandle = function (x, y) {
    for (var i = 0; i < this.handles.length; i++) {
      if (this.handles[i].contains(x, y)) return i;
    }

    return -1;
  };

  Board.prototype.drawHandle = function (handle, context) {
    context.strokeStyle = "yellow";
    context.beginPath();
    context.rect(handle.x, handle.y, handle.width, handle.height);
    context.stroke();
    context.closePath();
  };

  return Board;
}();

exports.Board = Board;
},{}],"components/Mcu.ts":[function(require,module,exports) {
"use strict";
/*
 * Microcontroller unit (MCU)
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Board_1 = require("./Board");

var Mcu =
/** @class */
function (_super) {
  __extends(Mcu, _super);

  function Mcu() {
    return _super !== null && _super.apply(this, arguments) || this;
  }

  return Mcu;
}(Board_1.Board);

exports.Mcu = Mcu;
},{"./Board":"components/Board.ts"}],"math/Rectangle.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Rectangle =
/** @class */
function () {
  function Rectangle(x, y, width, height) {
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  Rectangle.prototype.getXmax = function () {
    return this.x + this.width;
  };

  Rectangle.prototype.getYmax = function () {
    return this.y + this.height;
  };

  Rectangle.prototype.getXmin = function () {
    return this.x;
  };

  Rectangle.prototype.getYmin = function () {
    return this.y;
  };

  Rectangle.prototype.getCenterX = function () {
    return this.x + this.width / 2;
  };

  Rectangle.prototype.getCenterY = function () {
    return this.y + this.height / 2;
  };

  Rectangle.prototype.contains = function (px, py) {
    return px > this.x && px < this.x + this.width && py > this.y && py < this.y + this.height;
  };

  Rectangle.prototype.intersectRect = function (r) {
    return !(r.x > this.getXmax() || r.getXmax() < this.x || r.y > this.getYmax() || r.getYmax() < this.y);
  };

  Rectangle.prototype.toString = function () {
    return "[" + this.x + ", " + this.y + ", " + this.width + ", " + this.height + "]";
  };

  return Rectangle;
}();

exports.Rectangle = Rectangle;
},{}],"img/raspberry-pi.png":[function(require,module,exports) {
module.exports = "/raspberry-pi.73a9a4c7.png";
},{}],"components/RaspberryPi.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Mcu_1 = require("./Mcu");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // this is needed as we use augmented methods of canvas defined in Main.ts
// @ts-ignore


var raspberry_pi_png_1 = __importDefault(require("../img/raspberry-pi.png"));

var RaspberryPi =
/** @class */
function (_super) {
  __extends(RaspberryPi, _super);

  function RaspberryPi(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();
      var menu = document.getElementById(Main_1.contextMenus.raspberryPi.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      Main_1.contextMenus.raspberryPi.raspberryPi = _this;
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      // e.preventDefault();
      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(395, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(395, 245, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 245, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = raspberry_pi_png_1.default;

    _this.setX(20);

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  RaspberryPi.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 10;
    ctx.shadowOffsetY = 10;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 10;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  RaspberryPi.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  RaspberryPi.prototype.updateFirebase = function (value) {// TODO
  };

  RaspberryPi.prototype.updateFromFirebase = function () {// TODO
  };

  return RaspberryPi;
}(Mcu_1.Mcu);

exports.RaspberryPi = RaspberryPi;
},{"./Mcu":"components/Mcu.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/raspberry-pi.png":"img/raspberry-pi.png"}],"components/Hat.ts":[function(require,module,exports) {
"use strict";
/*
 * Hardware attached on top (HAT) for Raspberry Pi, equivalent to shields for Arduino
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Board_1 = require("./Board");

var RaspberryPi_1 = require("./RaspberryPi");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main");

var Hat =
/** @class */
function (_super) {
  __extends(Hat, _super);

  function Hat() {
    return _super !== null && _super.apply(this, arguments) || this;
  }

  Hat.prototype.attach = function (raspberryPi) {
    if (raspberryPi != null) {
      raspberryPi.hat = this;
    } else {
      if (this.raspberryPi != null) {
        this.raspberryPi.hat = null;
      }
    }

    this.raspberryPi = raspberryPi;

    if (raspberryPi != null) {
      this.setX(raspberryPi.getX());
      this.setY(raspberryPi.getY());
      localStorage.setItem("Attachment: " + this.getUid(), raspberryPi.uid);
    } else {
      localStorage.removeItem("Attachment: " + this.getUid());
    }
  };

  Hat.prototype.whichRaspberryPi = function () {
    var r1 = new Rectangle_1.Rectangle(this.getX(), this.getY(), this.getWidth(), this.getHeight());

    for (var i = Main_1.system.mcus.length - 1; i >= 0; i--) {
      var mcu = Main_1.system.mcus[i];

      if (mcu instanceof RaspberryPi_1.RaspberryPi) {
        var a = 20;

        if (r1.intersectRect(new Rectangle_1.Rectangle(mcu.getX() + a, mcu.getY() + a, mcu.getWidth() - 2 * a, mcu.getHeight() - 2 * a))) {
          return i;
        }
      }
    }

    return -1;
  };

  return Hat;
}(Board_1.Board);

exports.Hat = Hat;
},{"./Board":"components/Board.ts","./RaspberryPi":"components/RaspberryPi.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts"}],"components/LedDisplay.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var LedDisplay =
/** @class */
function () {
  function LedDisplay(board, name, x, y, width, height) {
    this.color = "lightgreen";
    this.fontSize = "70px";
    this.fontFamily = "digital-clock-font";
    this.board = board;
    this.name = name;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  LedDisplay.prototype.draw = function (ctx) {
    ctx.save();

    if (this.character) {
      ctx.fillStyle = "lightgreen";
      ctx.font = this.fontSize + " " + this.fontFamily;
      ctx.fillText(this.character, this.x, this.y);
    }

    ctx.restore();
  };

  LedDisplay.prototype.setCharacter = function (character) {
    if (character != null && character.length != 1) throw "Only one character is allowed";
    this.character = character;
  };

  LedDisplay.prototype.contains = function (x, y) {
    return x > this.x && x < this.x + this.width && y > this.y && y < this.y + this.height;
  };

  return LedDisplay;
}();

exports.LedDisplay = LedDisplay;
},{}],"components/Buzzer.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Buzzer =
/** @class */
function () {
  function Buzzer(board, name, x, y, width, height) {
    this.board = board;
    this.name = name;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  Buzzer.prototype.draw = function (ctx) {
    ctx.save();

    if (this.on) {
      ctx.lineWidth = 3;
      ctx.strokeStyle = "white";
      var centerX = this.x + this.width / 2;
      var centerY = this.y + this.height / 2;
      var angle = 0.25 * Math.PI;

      for (var i = 0; i < 4; i++) {
        var r = 5 + i * 5;
        ctx.beginPath();
        ctx.arc(centerX, centerY, r, -angle, angle);
        ctx.stroke();
        ctx.beginPath();
        ctx.arc(centerX, centerY, r, Math.PI - angle, Math.PI + angle);
        ctx.stroke();
      }
    }

    ctx.restore();
  };

  Buzzer.prototype.contains = function (x, y) {
    return x > this.x && x < this.x + this.width && y > this.y && y < this.y + this.height;
  };

  Buzzer.prototype.beepButton = function (button) {
    switch (button) {
      case "A":
        this.beep(1, 800, 200);
        break;

      case "B":
        this.beep(1, 400, 200);
        break;

      case "C":
        this.beep(1, 200, 200);
        break;
    }
  };

  Buzzer.prototype.beep = function (volume, frequency, duration) {
    this.on = true;
    if (this.audioContext == null) this.audioContext = new AudioContext();
    var v = this.audioContext.createOscillator();
    var u = this.audioContext.createGain();
    v.connect(u);
    v.frequency.value = frequency;
    v.type = "square";
    u.connect(this.audioContext.destination);
    u.gain.value = volume * 0.01;
    v.start(this.audioContext.currentTime);
    v.stop(this.audioContext.currentTime + duration * 0.001);
    var that = this;
    setTimeout(function () {
      that.on = false;
      that.board.draw();
    }, 200);
  };

  return Buzzer;
}();

exports.Buzzer = Buzzer;
},{}],"Util.ts":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
/*
 * @author Charles Xie
 */

var Util =
/** @class */
function () {
  function Util() {}

  Util.countDigits = function (x) {
    return x.toString().length;
  };

  Util.adjust = function (color, amount) {
    return '#' + color.replace(/^#/, '').replace(/../g, function (color) {
      return ('0' + Math.min(255, Math.max(0, parseInt(color, 16) + amount)).toString(16)).substr(-2);
    });
  };

  Util.hexToRgb = function (hex) {
    var result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
    return result ? {
      r: parseInt(result[1], 16),
      g: parseInt(result[2], 16),
      b: parseInt(result[3], 16)
    } : null;
  };

  Util.rgbToHex = function (r, g, b) {
    return "#" + ((1 << 24) + (r << 16) + (g << 8) + b).toString(16).slice(1);
  };

  Util.rgbToHue = function (r, g, b) {
    r /= 255;
    g /= 255;
    b /= 255;
    var max = Math.max(r, g, b);
    var min = Math.min(r, g, b);
    var c = max - min;
    var hue;

    if (c == 0) {
      hue = 0;
    } else {
      switch (max) {
        case r:
          var segment = (g - b) / c;
          var shift = 0 / 60; // R° / (360° / hex sides)

          if (segment < 0) {
            // hue > 180, full rotation
            shift = 360 / 60; // R° / (360° / hex sides)
          }

          hue = segment + shift;
          break;

        case g:
          var segment = (b - r) / c;
          var shift = 120 / 60; // G° / (360° / hex sides)

          hue = segment + shift;
          break;

        case b:
          var segment = (r - g) / c;
          var shift = 240 / 60; // B° / (360° / hex sides)

          hue = segment + shift;
          break;
      }
    }

    return hue * 60; // hue is in [0,6], scale it up
  };

  Util.rgbToHsl = function (r, g, b) {
    var r1 = r / 255;
    var g1 = g / 255;
    var b1 = b / 255;
    var maxColor = Math.max(r1, g1, b1);
    var minColor = Math.min(r1, g1, b1);
    var lightness = (maxColor + minColor) / 2;
    var saturation = 0;
    var hue = 0;

    if (maxColor != minColor) {
      if (lightness < 0.5) {
        saturation = (maxColor - minColor) / (maxColor + minColor);
      } else {
        saturation = (maxColor - minColor) / (2 - maxColor - minColor);
      }

      if (r1 == maxColor) {
        hue = (g1 - b1) / (maxColor - minColor);
      } else if (g1 == maxColor) {
        hue = 2.0 + (b1 - r1) / (maxColor - minColor);
      } else {
        hue = 4.0 + (r1 - g1) / (maxColor - minColor);
      }
    }

    lightness = lightness * 100;
    saturation = saturation * 100;
    hue = hue * 60;

    if (hue < 0) {
      hue += 360;
    }

    return {
      h: hue,
      s: saturation,
      l: lightness
    };
  };
  /*
   * Converts an HSL color value to RGB. Conversion formula
   * adapted from http://en.wikipedia.org/wiki/HSL_color_space.
   * Assumes h, s, and l are contained in the set [0, 1] and
   * returns r, g, and b in the set [0, 255].
   *
   * @param   {number}  h       The hue
   * @param   {number}  s       The saturation
   * @param   {number}  l       The lightness
   * @return  {Array}           The RGB representation
   */


  Util.hslToRgb = function (h, s, l) {
    var r, g, b;

    if (s == 0) {
      r = g = b = l; // achromatic
    } else {
      var hue2rgb = function hue2rgb(p, q, t) {
        if (t < 0) t += 1;
        if (t > 1) t -= 1;
        if (t < 1 / 6) return p + (q - p) * 6 * t;
        if (t < 1 / 2) return q;
        if (t < 2 / 3) return p + (q - p) * (2 / 3 - t) * 6;
        return p;
      };

      var q = l < 0.5 ? l * (1 + s) : l + s - l * s;
      var p = 2 * l - q;
      r = hue2rgb(p, q, h + 1 / 3);
      g = hue2rgb(p, q, h);
      b = hue2rgb(p, q, h - 1 / 3);
    }

    return {
      r: Math.round(r * 255),
      g: Math.round(g * 255),
      b: Math.round(b * 255)
    };
  };

  Util.getHueColor = function (color) {
    var c = Util.hexToRgb(color);
    c = Util.hslToRgb(Util.rgbToHue(c.r, c.g, c.b) / 360, 0.5, 0.5);
    return Util.rgbToHex(c.r, c.g, c.b);
  }; // expects 0 <= h, s, v <= 1


  Util.hsvToRgb = function (h, s, v) {
    var r, g, b;
    var i = Math.floor(h * 6);
    var f = h * 6 - i;
    var p = v * (1 - s);
    var q = v * (1 - f * s);
    var t = v * (1 - (1 - f) * s);

    switch (i % 6) {
      case 0:
        r = v, g = t, b = p;
        break;

      case 1:
        r = q, g = v, b = p;
        break;

      case 2:
        r = p, g = v, b = t;
        break;

      case 3:
        r = p, g = q, b = v;
        break;

      case 4:
        r = t, g = p, b = v;
        break;

      case 5:
        r = v, g = p, b = q;
        break;
    }

    return {
      r: Math.round(r * 255),
      g: Math.round(g * 255),
      b: Math.round(b * 255)
    };
  }; // output 0 <= h, s, v <= 1


  Util.rgbToHsv = function (r, g, b) {
    var max = Math.max(r, g, b);
    var min = Math.min(r, g, b);
    var d = max - min;
    var h;
    var s = max === 0 ? 0 : d / max;
    var v = max / 255;

    switch (max) {
      case min:
        h = 0;
        break;

      case r:
        h = g - b + d * (g < b ? 6 : 0);
        h /= 6 * d;
        break;

      case g:
        h = b - r + d * 2;
        h /= 6 * d;
        break;

      case b:
        h = r - g + d * 4;
        h /= 6 * d;
        break;
    }

    return {
      h: h,
      s: s,
      v: v
    };
  };

  return Util;
}();

exports.Util = Util;
},{}],"components/LedLight.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Util_1 = require("../Util");

var LedLight =
/** @class */
function () {
  function LedLight(board, name, color, radius, rays, rayLength, x, y, width, height) {
    this.color = "#ff0000";
    this.radius = 4;
    this.rays = 8;
    this.rayLength = 5;
    this.on = false;
    this.board = board;
    this.name = name;
    this.color = color;
    this.radius = radius;
    this.rays = rays;
    this.rayLength = rayLength;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  LedLight.prototype.draw = function (ctx) {
    ctx.save();

    if (this.on) {
      ctx.beginPath();
      ctx.lineWidth = 1;
      var centerX = this.x + this.width / 2;
      var centerY = this.y + this.height / 2;

      if (this.radius < 8) {
        var x1 = void 0,
            y1 = void 0,
            x2 = void 0,
            y2 = void 0;
        var angle = void 0,
            cos = void 0,
            sin = void 0;
        ctx.strokeStyle = "white";

        for (var i = 0; i < this.rays; i++) {
          angle = i * Math.PI * 2 / this.rays;
          cos = Math.cos(angle);
          sin = Math.sin(angle);
          x1 = centerX + this.radius * cos * 2 / 3;
          y1 = centerY + this.radius * sin * 2 / 3;
          x2 = centerX + (this.radius + this.rayLength) * cos;
          y2 = centerY + (this.radius + this.rayLength) * sin;
          ctx.beginPath();
          ctx.moveTo(x1, y1);
          ctx.lineTo(x2, y2);
          ctx.stroke();
          ctx.closePath();
        }

        ctx.fillStyle = this.color;
      } else {
        var gradient = ctx.createRadialGradient(centerX, centerY, 0, centerX, centerY, this.radius);
        gradient.addColorStop(1, "rgba(255, 255, 255, 0)");
        gradient.addColorStop(0.5, Util_1.Util.adjust(this.color, -100));
        gradient.addColorStop(0, Util_1.Util.adjust(this.color, 100));
        ctx.fillStyle = gradient;
      }

      ctx.arc(centerX, centerY, this.radius, 0, Math.PI * 2);
      ctx.closePath();
      ctx.fill();
    } else {
      ctx.beginPath();
      ctx.lineWidth = 2;
      ctx.strokeStyle = this.color;
      ctx.arc(this.x + this.width / 2, this.y + this.height / 2, this.radius, 0, Math.PI * 2);
      ctx.stroke();
      ctx.closePath();
    }

    ctx.restore();
  };

  LedLight.prototype.contains = function (x, y) {
    return x > this.x && x < this.x + this.width && y > this.y && y < this.y + this.height;
  };

  LedLight.prototype.toggle = function (x, y) {
    var inside = this.contains(x, y);

    if (inside) {
      this.on = !this.on;
    }

    return inside;
  };

  return LedLight;
}();

exports.LedLight = LedLight;
},{"../Util":"Util.ts"}],"components/Button.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Button =
/** @class */
function () {
  function Button(board, name, x, y, width, height) {
    this.on = false;
    this.pressedColor = '#66cccccc';
    this.board = board;
    this.name = name;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  Button.prototype.draw = function (ctx) {
    if (this.on) {
      ctx.fillStyle = this.pressedColor;
      ctx.fillRect(this.x, this.y, this.width, this.height);
    }

    ctx.lineWidth = 1;
    ctx.strokeStyle = 'black';
    ctx.rect(this.x, this.y, this.width, this.height);
    ctx.stroke();
  };

  Button.prototype.contains = function (x, y) {
    return x > this.x && x < this.x + this.width && y > this.y && y < this.y + this.height;
  };

  return Button;
}();

exports.Button = Button;
},{}],"components/Sensor.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Sensor =
/** @class */
function () {
  function Sensor(board, name, unit, x, y, width, height) {
    this.data = [];
    this.collectionInterval = 1; // in seconds

    this.on = false;
    this.pressedColor = 'white';
    this.board = board;
    this.name = name;
    this.unit = unit;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
  }

  Sensor.prototype.draw = function (ctx) {
    if (this.on) {
      ctx.fillStyle = this.pressedColor;
      ctx.fillRect(this.x, this.y, this.width, this.height);
    }

    ctx.lineWidth = 1;
    ctx.strokeStyle = 'black';
    ctx.rect(this.x, this.y, this.width, this.height);
    ctx.stroke();
  };

  Sensor.prototype.contains = function (x, y) {
    return x > this.x && x < this.x + this.width && y > this.y && y < this.y + this.height;
  };

  return Sensor;
}();

exports.Sensor = Sensor;
},{}],"tools/ColorPicker.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Util_1 = require("../Util");

var ColorPicker =
/** @class */
function () {
  function ColorPicker() {
    this.rgbaColor = "#ffffffff";
    this.savedBlockX = -10;
    this.savedBlockY = -10;
    this.savedStripY = -10;
    this.colorBlock = document.getElementById('color-block');
    this.ctx1 = this.colorBlock.getContext('2d');
    this.colorStrip = document.getElementById('color-strip');
    this.ctx2 = this.colorStrip.getContext('2d');
    this.colorBlock.addEventListener("mousedown", this.mouseDownBlock.bind(this), false);
    this.colorBlock.addEventListener("mouseup", this.mouseUpBlock.bind(this), false);
    this.colorStrip.addEventListener("click", this.mouseClickStrip.bind(this), false);
    this.colorStrip.addEventListener("mousedown", this.mouseDownStrip.bind(this), false);
    this.colorStrip.addEventListener("mouseup", this.mouseUpStrip.bind(this), false);
  }

  ColorPicker.prototype.setColorLabel = function (colorLabel) {
    this.colorLabel = colorLabel;
  };

  ColorPicker.prototype.setColorCode = function (colorCode) {
    this.colorCode = colorCode;
  };

  ColorPicker.prototype.setSelectedColor = function (color) {
    this.rgbaColor = color;
    var c = Util_1.Util.hexToRgb(color);
    this.savedStripY = this.colorStrip.height * Util_1.Util.rgbToHue(c.r, c.g, c.b) / 360;

    if (this.colorLabel) {
      this.colorLabel.style.backgroundColor = color;
    }

    if (this.colorCode) {
      this.colorCode.value = color;
      this.colorCode.select();
    }
  }; // TODO: the position should be calculated by inverting the dual-axis shading used in fillGradient()
  // right now, we just use a brute-force method to search for the point that maps to the selected color


  ColorPicker.prototype.setSelectedPoint = function () {
    var c = Util_1.Util.hexToRgb(this.rgbaColor);
    var imageData = this.ctx1.getImageData(0, 0, this.colorBlock.width, this.colorBlock.height).data;
    var n = imageData.length / 4;

    for (var i = 0; i < n; i++) {
      if (imageData[4 * i] == c.r && imageData[4 * i + 1] == c.g && imageData[4 * i + 2] == c.b) {
        this.savedBlockX = i % this.colorBlock.width;
        this.savedBlockY = i / this.colorBlock.width;
        this.draw();
        break;
      }
    }
  };

  ColorPicker.prototype.getSelectedColor = function () {
    return this.rgbaColor;
  };

  ColorPicker.prototype.draw = function () {
    // draw color block
    this.ctx1.clearRect(0, 0, this.colorBlock.width, this.colorBlock.height);
    this.ctx1.rect(0, 0, this.colorBlock.width, this.colorBlock.height);
    this.fillGradient();
    this.ctx1.strokeStyle = "white";
    this.ctx1.lineWidth = 2;
    this.ctx1.beginPath();
    this.ctx1.arc(this.savedBlockX, this.savedBlockY, 6, 0, 2 * Math.PI);
    this.ctx1.stroke(); // draw color strip

    this.ctx2.clearRect(0, 0, this.colorStrip.width, this.colorStrip.height);
    this.ctx2.rect(0, 0, this.colorStrip.width, this.colorStrip.height);
    var gradient = this.ctx2.createLinearGradient(0, 0, 0, this.colorStrip.height);
    gradient.addColorStop(0, '#ff0000');
    gradient.addColorStop(1.0 / 6.0, '#ffff00');
    gradient.addColorStop(2.0 / 6.0, '#00ff00');
    gradient.addColorStop(3.0 / 6.0, '#00ffff');
    gradient.addColorStop(4.0 / 6.0, '#0000ff');
    gradient.addColorStop(5.0 / 6.0, '#ff00ff');
    gradient.addColorStop(1, '#ff0000');
    this.ctx2.fillStyle = gradient;
    this.ctx2.fill();
    this.ctx2.fillStyle = "white";
    this.ctx2.beginPath();
    this.ctx2.rect(0, this.savedStripY, this.colorStrip.width, 6);
    this.ctx2.fill();
    this.ctx2.strokeStyle = "black";
    this.ctx2.lineWidth = 2;
    this.ctx2.beginPath();
    this.ctx2.rect(0, this.savedStripY, this.colorStrip.width, 8);
    this.ctx2.stroke();
  };

  ColorPicker.prototype.fillGradient = function () {
    this.ctx1.fillStyle = Util_1.Util.getHueColor(this.rgbaColor);
    this.ctx1.fillRect(0, 0, this.colorBlock.width, this.colorBlock.height);
    var grdWhite = this.ctx1.createLinearGradient(0, 0, this.colorBlock.width, 0);
    grdWhite.addColorStop(0, 'rgba(255, 255, 255, 1)');
    grdWhite.addColorStop(1, 'rgba(255, 255, 255, 0)');
    this.ctx1.fillStyle = grdWhite;
    this.ctx1.fillRect(0, 0, this.colorBlock.width, this.colorBlock.height);
    var grdBlack = this.ctx1.createLinearGradient(0, 0, 0, this.colorBlock.height);
    grdBlack.addColorStop(0, 'rgba(0, 0, 0, 0)');
    grdBlack.addColorStop(1, 'rgba(0, 0, 0, 1)');
    this.ctx1.fillStyle = grdBlack;
    this.ctx1.fillRect(0, 0, this.colorBlock.width, this.colorBlock.height);
  };

  ColorPicker.prototype.mouseDownStrip = function (e) {
    e.stopPropagation();
  };

  ColorPicker.prototype.mouseUpStrip = function (e) {
    e.stopPropagation();
  };

  ColorPicker.prototype.mouseClickStrip = function (e) {
    this.savedStripY = e.offsetY;
    var imageData = this.ctx2.getImageData(e.offsetX, e.offsetY, 1, 1).data;
    this.rgbaColor = Util_1.Util.rgbToHex(imageData[0], imageData[1], imageData[2]);
    this.draw();
    e.stopPropagation();
  };

  ColorPicker.prototype.mouseDownBlock = function (e) {
    this.changeColor(e);
    e.stopPropagation();
  };

  ColorPicker.prototype.mouseUpBlock = function (e) {
    this.savedBlockX = e.offsetX;
    this.savedBlockY = e.offsetY;
    this.draw();
    e.stopPropagation();
  };

  ColorPicker.prototype.changeColor = function (e) {
    var imageData = this.ctx1.getImageData(e.offsetX, e.offsetY, 1, 1).data;
    this.rgbaColor = Util_1.Util.rgbToHex(imageData[0], imageData[1], imageData[2]);

    if (this.colorLabel) {
      this.colorLabel.style.backgroundColor = this.rgbaColor;
    }

    if (this.colorCode) {
      this.colorCode.value = this.rgbaColor;
    }
  };

  return ColorPicker;
}();

exports.ColorPicker = ColorPicker;
},{"../Util":"Util.ts"}],"img/rainbow-hat.png":[function(require,module,exports) {
module.exports = "/rainbow-hat.151dd447.png";
},{}],"components/RainbowHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Rainbow HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var LedDisplay_1 = require("./LedDisplay");

var Buzzer_1 = require("./Buzzer");

var LedLight_1 = require("./LedLight");

var Button_1 = require("./Button");

var Sensor_1 = require("./Sensor");

var System_1 = require("../System");

var Util_1 = require("../Util");

var Rectangle_1 = require("../math/Rectangle");

var ColorPicker_1 = require("../tools/ColorPicker");

var Main_1 = require("../Main"); // @ts-ignore


var rainbow_hat_png_1 = __importDefault(require("../img/rainbow-hat.png"));

var RainbowHat =
/** @class */
function (_super) {
  __extends(RainbowHat, _super);

  function RainbowHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.rgbLedLights = [];
    _this.alphanumericDisplays = [];
    _this.decimalPointDisplays = [];
    _this.indexOfSelectedRgbLedLight = -1;
    _this.stateId = "rainbow_hat_default";

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      _this.indexOfSelectedRgbLedLight = -1;

      for (var i = 0; i < _this.rgbLedLights.length; i++) {
        if (_this.rgbLedLights[i].contains(dx, dy)) {
          _this.indexOfSelectedRgbLedLight = i;
          break;
        }
      }

      if (_this.indexOfSelectedRgbLedLight >= 0) {
        Main_1.contextMenus.colorPicker.rainbowHat = _this;
        var menu = document.getElementById(Main_1.contextMenus.colorPicker.id);
        menu.style.left = e.clientX + "px";
        menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
        menu.classList.add("show-menu");

        if (Main_1.system.colorPicker == null) {
          Main_1.system.colorPicker = new ColorPicker_1.ColorPicker();
        }

        Main_1.system.colorPicker.setColorLabel(document.getElementById("colorpicker-label"));
        Main_1.system.colorPicker.setColorCode(document.getElementById("colorpicker-hex-code"));
        Main_1.system.colorPicker.setSelectedColor(_this.rgbLedLights[_this.indexOfSelectedRgbLedLight].color);
        Main_1.system.colorPicker.draw();
        Main_1.system.colorPicker.setSelectedPoint();
        document.getElementById("colorpicker-title").innerText = "RGB LED Light " + _this.indexOfSelectedRgbLedLight;
      } else {
        Main_1.contextMenus.rainbowHat.hat = _this;
        var menu = document.getElementById(Main_1.contextMenus.rainbowHat.id);
        menu.style.left = e.clientX + "px";
        menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
        menu.classList.add("show-menu");
        var attachMenuItem = document.getElementById("rainbow-hat-attach-menu-item");
        var detachMenuItem = document.getElementById("rainbow-hat-detach-menu-item");

        if (_this.raspberryPi != null) {
          attachMenuItem.className = "menu-item disabled";
          detachMenuItem.className = "menu-item";
        } else {
          var i = _this.whichRaspberryPi();

          attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
          detachMenuItem.className = "menu-item disabled";
        }
      }
    };

    _this.mouseDown = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      var context = _this.canvas.getContext("2d");

      if (_this.redLedLight.toggle(dx, dy)) {
        _this.updateFirebase({
          redLed: _this.redLedLight.on
        });

        _this.redLedLight.draw(context);

        return;
      }

      if (_this.greenLedLight.toggle(dx, dy)) {
        _this.updateFirebase({
          greenLed: _this.greenLedLight.on
        });

        _this.greenLedLight.draw(context);

        return;
      }

      if (_this.blueLedLight.toggle(dx, dy)) {
        _this.updateFirebase({
          blueLed: _this.blueLedLight.on
        });

        _this.blueLedLight.draw(context);

        return;
      }

      if (_this.buttonA.contains(dx, dy)) {
        _this.buttonA.on = true;

        _this.buttonA.draw(context);

        _this.redLedLight.on = true;

        _this.redLedLight.draw(context);

        _this.updateFirebase({
          redLed: true
        });

        _this.buzzer.beepButton("A");

        return;
      }

      if (_this.buttonB.contains(dx, dy)) {
        _this.buttonB.on = true;

        _this.buttonB.draw(context);

        _this.greenLedLight.on = true;

        _this.greenLedLight.draw(context);

        _this.updateFirebase({
          greenLed: true
        });

        _this.buzzer.beepButton("B");

        return;
      }

      if (_this.buttonC.contains(dx, dy)) {
        _this.buttonC.on = true;

        _this.buttonC.draw(context);

        _this.blueLedLight.on = true;

        _this.blueLedLight.draw(context);

        _this.updateFirebase({
          blueLed: true
        });

        _this.buzzer.beepButton("C");

        return;
      }
    };

    _this.mouseUp = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      var context = _this.canvas.getContext("2d");

      if (_this.buttonA.contains(dx, dy)) {
        _this.buttonA.on = false;

        _this.buttonA.draw(context);

        _this.redLedLight.on = false;

        _this.redLedLight.draw(context);

        _this.updateFirebase({
          redLed: false
        });

        return;
      }

      if (_this.buttonB.contains(dx, dy)) {
        _this.buttonB.on = false;

        _this.buttonB.draw(context);

        _this.greenLedLight.on = false;

        _this.greenLedLight.draw(context);

        _this.updateFirebase({
          greenLed: false
        });

        return;
      }

      if (_this.buttonC.contains(dx, dy)) {
        _this.buttonC.on = false;

        _this.buttonC.draw(context);

        _this.blueLedLight.on = false;

        _this.blueLedLight.draw(context);

        _this.updateFirebase({
          blueLed: false
        });

        return;
      }

      if (_this.temperatureSensor.contains(dx, dy)) {
        Main_1.system.temperatureGraph.setVisible(!Main_1.system.temperatureGraph.isVisible());

        if (Main_1.system.temperatureGraph.isVisible()) {
          Main_1.system.temperatureGraph.draw();
          Main_1.system.temperatureGraph.bringForward();
        }

        localStorage.setItem("Visible: " + Main_1.system.temperatureGraph.getUid(), Main_1.system.temperatureGraph.isVisible() ? "true" : "false");
        return;
      }

      if (_this.barometricPressureSensor.contains(dx, dy)) {
        Main_1.system.pressureGraph.setVisible(!Main_1.system.pressureGraph.isVisible());

        if (Main_1.system.pressureGraph.isVisible()) {
          Main_1.system.pressureGraph.draw();
          Main_1.system.pressureGraph.bringForward();
        }

        localStorage.setItem("Visible: " + Main_1.system.pressureGraph.getUid(), Main_1.system.pressureGraph.isVisible() ? "true" : "false");
        return;
      }
    };

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.redLedLight.contains(dx, dy)) {
        _this.mouseOverObject = _this.redLedLight;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.greenLedLight.contains(dx, dy)) {
        _this.mouseOverObject = _this.greenLedLight;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.blueLedLight.contains(dx, dy)) {
        _this.mouseOverObject = _this.blueLedLight;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.buttonA.contains(dx, dy)) {
        _this.mouseOverObject = _this.buttonA;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.buttonB.contains(dx, dy)) {
        _this.mouseOverObject = _this.buttonB;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.buttonC.contains(dx, dy)) {
        _this.mouseOverObject = _this.buttonC;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.temperatureSensor.contains(dx, dy)) {
        _this.mouseOverObject = _this.temperatureSensor;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.barometricPressureSensor.contains(dx, dy)) {
        _this.mouseOverObject = _this.barometricPressureSensor;
        _this.canvas.style.cursor = "pointer";
      } else if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        var onRgbLedLight = false;

        for (var i = 0; i < _this.rgbLedLights.length; i++) {
          if (_this.rgbLedLights[i].contains(dx, dy)) {
            _this.mouseOverObject = _this.rgbLedLights[i];
            _this.canvas.style.cursor = "pointer";
            onRgbLedLight = true;
            break;
          }
        }

        if (!onRgbLedLight) {
          _this.mouseOverObject = null;
          _this.canvas.style.cursor = "default";
        }
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.redLedLight = new LedLight_1.LedLight(_this, "LED Light", "red", 4, 8, 10, 65, 233, 18, 8);
    _this.greenLedLight = new LedLight_1.LedLight(_this, "LED Light", "green", 4, 8, 10, 147, 233, 18, 8);
    _this.blueLedLight = new LedLight_1.LedLight(_this, "LED Light", "blue", 4, 8, 10, 230, 233, 18, 8);
    _this.buttonA = new Button_1.Button(_this, "Button A", 38, 245, 72, 24);
    _this.buttonB = new Button_1.Button(_this, "Button B", 120, 245, 72, 24);
    _this.buttonC = new Button_1.Button(_this, "Button C", 203, 245, 72, 24);
    _this.buzzer = new Buzzer_1.Buzzer(_this, "Piezo Buzzer", 35, 170, 20, 20);
    _this.temperatureSensor = new Sensor_1.Sensor(_this, "Temperature", "°C", 152, 108, 12, 12);
    _this.barometricPressureSensor = new Sensor_1.Sensor(_this, "Pressure", "hPa", 187, 115, 20, 10);

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 251, 78, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 218, 62, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 183, 53, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 147, 50, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 111, 53, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 76, 62, 20, 20));

    _this.rgbLedLights.push(new LedLight_1.LedLight(_this, "RGB LED Light", "black", 16, 12, 2, 44, 78, 20, 20));

    _this.alphanumericDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 34, 214, 33, 65));

    _this.alphanumericDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 95, 214, 33, 65));

    _this.alphanumericDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 156, 214, 33, 65));

    _this.alphanumericDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 218, 214, 33, 65));

    _this.decimalPointDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 78, 214, 33, 65));

    _this.decimalPointDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 139, 214, 33, 65));

    _this.decimalPointDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 200, 214, 33, 65));

    _this.decimalPointDisplays.push(new LedDisplay_1.LedDisplay(_this, "LED Display", 261, 214, 33, 65));

    for (var i = 0; i < _this.decimalPointDisplays.length; i++) {
      _this.decimalPointDisplays[i].fontSize = "40px";
    }

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = rainbow_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  RainbowHat.prototype.setSelectedRgbLedLightColor = function (color) {
    if (this.indexOfSelectedRgbLedLight >= 0) {
      this.rgbLedLights[this.indexOfSelectedRgbLedLight].color = color;
      this.draw();
      var list = [];

      for (var i = 0; i < this.rgbLedLights.length; i++) {
        var a = [];
        var c = Util_1.Util.hexToRgb(this.rgbLedLights[i].color);
        a.push(c.r);
        a.push(c.g);
        a.push(c.b);
        list.push(a);
      }

      this.updateFirebase({
        rainbowRgb: list
      });
    }
  };

  RainbowHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.redLedLight.draw(ctx);
    this.greenLedLight.draw(ctx);
    this.blueLedLight.draw(ctx);
    this.buttonA.draw(ctx);
    this.buttonB.draw(ctx);
    this.buttonC.draw(ctx);
    this.buzzer.draw(ctx);
    this.temperatureSensor.draw(ctx);
    this.barometricPressureSensor.draw(ctx);
    this.drawToolTips();

    for (var i = 0; i < this.rgbLedLights.length; i++) {
      this.rgbLedLights[i].draw(ctx);
    }

    for (var i = 0; i < this.alphanumericDisplays.length; i++) {
      this.alphanumericDisplays[i].draw(ctx);
    }

    for (var i = 0; i < this.decimalPointDisplays.length; i++) {
      this.decimalPointDisplays[i].draw(ctx);
    }
  };

  RainbowHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.redLedLight) {
      x += this.redLedLight.x + this.redLedLight.width / 2;
      y += this.redLedLight.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Red LED light', true);
    } else if (this.mouseOverObject == this.greenLedLight) {
      x += this.greenLedLight.x + this.greenLedLight.width / 2;
      y += this.greenLedLight.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Green LED light', true);
    } else if (this.mouseOverObject == this.blueLedLight) {
      x += this.blueLedLight.x + this.blueLedLight.width / 2;
      y += this.blueLedLight.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Blue LED light', true);
    } else if (this.mouseOverObject == this.buttonA) {
      x += this.buttonA.x + this.buttonA.width / 2;
      y += this.buttonA.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Button A', true);
    } else if (this.mouseOverObject == this.buttonB) {
      x += this.buttonB.x + this.buttonB.width / 2;
      y += this.buttonB.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Button B', true);
    } else if (this.mouseOverObject == this.buttonC) {
      x += this.buttonC.x + this.buttonC.width / 2;
      y += this.buttonC.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Button C', true);
    } else if (this.mouseOverObject == this.temperatureSensor) {
      x += this.temperatureSensor.x + this.temperatureSensor.width / 2;
      y += this.temperatureSensor.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Temperature sensor', true);
    } else if (this.mouseOverObject == this.barometricPressureSensor) {
      x += this.barometricPressureSensor.x + this.barometricPressureSensor.width / 2;
      y += this.barometricPressureSensor.y;
      context.drawTooltip(x, y, 20, 8, 10, 'Barometric pressure sensor', true);
    } else if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  RainbowHat.prototype.updateFirebase = function (value) {
    System_1.System.database.ref(this.stateId).update(value);
  }; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  RainbowHat.prototype.updateFromFirebase = function () {
    var that = this;
    System_1.System.database.ref().on("value", function (snapshot) {
      snapshot.forEach(function (child) {
        var childData = child.val();
        that.redLedLight.on = childData.redLed;
        that.greenLedLight.on = childData.greenLed;
        that.blueLedLight.on = childData.blueLed;

        if (that.redLedLight.on) {
          that.buzzer.beepButton("A");
        }

        if (that.greenLedLight.on) {
          that.buzzer.beepButton("B");
        }

        if (that.blueLedLight.on) {
          that.buzzer.beepButton("C");
        }

        if (childData.rainbowRgb) {
          for (var i = 0; i < that.rgbLedLights.length; i++) {
            var r = childData.rainbowRgb[i][0];
            var g = childData.rainbowRgb[i][1];
            var b = childData.rainbowRgb[i][2];
            that.rgbLedLights[i].on = r > 0 || g > 0 || b > 0;
            that.rgbLedLights[i].color = Util_1.Util.rgbToHex(r, g, b);
          }
        }

        if (childData.allowTemperatureTransmission) {
          that.temperatureSensor.collectionInterval = childData.sensorDataCollectionInterval ? childData.sensorDataCollectionInterval * 0.001 : 1;
          that.temperatureSensor.data.push(childData.temperature);
          Main_1.system.temperatureGraph.draw();
          var t = childData.temperature;
          var s = t.toString();
          var i = s.indexOf(".");

          if (i > 0 && i < 4) {
            that.decimalPointDisplays[i - 1].setCharacter(".");
          }

          var integerPart = s.substring(0, i);
          var decimalPart = s.substring(i + 1);
          s = (integerPart + decimalPart).substr(0, 4);

          for (i = 0; i < 4; i++) {
            that.alphanumericDisplays[i].setCharacter(s[i]);
          }
        }

        if (childData.allowBarometricPressureTransmission) {
          that.barometricPressureSensor.collectionInterval = childData.sensorDataCollectionInterval ? childData.sensorDataCollectionInterval * 0.001 : 1;
          that.barometricPressureSensor.data.push(childData.barometricPressure);
          Main_1.system.pressureGraph.draw();
          var t = childData.barometricPressure;
          var s = t.toString();
          var i = s.indexOf(".");

          if (i > 0 && i < 4) {
            that.decimalPointDisplays[i - 1].setCharacter(".");
          }

          var integerPart = s.substring(0, i);
          var decimalPart = s.substring(i + 1);
          s = (integerPart + decimalPart).substr(0, 4);

          for (i = 0; i < 4; i++) {
            that.alphanumericDisplays[i].setCharacter(s[i]);
          }
        }

        that.draw();
      });
    });
  };

  return RainbowHat;
}(Hat_1.Hat);

exports.RainbowHat = RainbowHat;
},{"./Hat":"components/Hat.ts","./LedDisplay":"components/LedDisplay.ts","./Buzzer":"components/Buzzer.ts","./LedLight":"components/LedLight.ts","./Button":"components/Button.ts","./Sensor":"components/Sensor.ts","../System":"System.ts","../Util":"Util.ts","../math/Rectangle":"math/Rectangle.ts","../tools/ColorPicker":"tools/ColorPicker.ts","../Main":"Main.ts","../img/rainbow-hat.png":"img/rainbow-hat.png"}],"img/sense-hat.png":[function(require,module,exports) {
module.exports = "/sense-hat.432da3f3.png";
},{}],"components/SenseHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Sense HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // @ts-ignore


var sense_hat_png_1 = __importDefault(require("../img/sense-hat.png"));

var SenseHat =
/** @class */
function (_super) {
  __extends(SenseHat, _super);

  function SenseHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      Main_1.contextMenus.senseHat.hat = _this;
      var menu = document.getElementById(Main_1.contextMenus.senseHat.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      var attachMenuItem = document.getElementById("sense-hat-attach-menu-item");
      var detachMenuItem = document.getElementById("sense-hat-detach-menu-item");

      if (_this.raspberryPi != null) {
        attachMenuItem.className = "menu-item disabled";
        detachMenuItem.className = "menu-item";
      } else {
        var i = _this.whichRaspberryPi();

        attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
        detachMenuItem.className = "menu-item disabled";
      }
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = sense_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  SenseHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  SenseHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  SenseHat.prototype.updateFirebase = function (value) {}; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  SenseHat.prototype.updateFromFirebase = function () {};

  return SenseHat;
}(Hat_1.Hat);

exports.SenseHat = SenseHat;
},{"./Hat":"components/Hat.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/sense-hat.png":"img/sense-hat.png"}],"img/capacitive-touch-hat.png":[function(require,module,exports) {
module.exports = "/capacitive-touch-hat.f7ea211b.png";
},{}],"components/CapacitiveTouchHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Capacitive Touch HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // @ts-ignore


var capacitive_touch_hat_png_1 = __importDefault(require("../img/capacitive-touch-hat.png"));

var CapacitiveTouchHat =
/** @class */
function (_super) {
  __extends(CapacitiveTouchHat, _super);

  function CapacitiveTouchHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      Main_1.contextMenus.capacitiveTouchHat.hat = _this;
      var menu = document.getElementById(Main_1.contextMenus.capacitiveTouchHat.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      var attachMenuItem = document.getElementById("capacitive-touch-hat-attach-menu-item");
      var detachMenuItem = document.getElementById("capacitive-touch-hat-detach-menu-item");

      if (_this.raspberryPi != null) {
        attachMenuItem.className = "menu-item disabled";
        detachMenuItem.className = "menu-item";
      } else {
        var i = _this.whichRaspberryPi();

        attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
        detachMenuItem.className = "menu-item disabled";
      }
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = capacitive_touch_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  CapacitiveTouchHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  CapacitiveTouchHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  CapacitiveTouchHat.prototype.updateFirebase = function (value) {}; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  CapacitiveTouchHat.prototype.updateFromFirebase = function () {};

  return CapacitiveTouchHat;
}(Hat_1.Hat);

exports.CapacitiveTouchHat = CapacitiveTouchHat;
},{"./Hat":"components/Hat.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/capacitive-touch-hat.png":"img/capacitive-touch-hat.png"}],"img/unicorn-hat.png":[function(require,module,exports) {
module.exports = "/unicorn-hat.9a311ef5.png";
},{}],"components/UnicornHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Unicorn HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // @ts-ignore


var unicorn_hat_png_1 = __importDefault(require("../img/unicorn-hat.png"));

var UnicornHat =
/** @class */
function (_super) {
  __extends(UnicornHat, _super);

  function UnicornHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      Main_1.contextMenus.unicornHat.hat = _this;
      var menu = document.getElementById(Main_1.contextMenus.unicornHat.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      var attachMenuItem = document.getElementById("unicorn-hat-attach-menu-item");
      var detachMenuItem = document.getElementById("unicorn-hat-detach-menu-item");

      if (_this.raspberryPi != null) {
        attachMenuItem.className = "menu-item disabled";
        detachMenuItem.className = "menu-item";
      } else {
        var i = _this.whichRaspberryPi();

        attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
        detachMenuItem.className = "menu-item disabled";
      }
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = unicorn_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  UnicornHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  UnicornHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  UnicornHat.prototype.updateFirebase = function (value) {}; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  UnicornHat.prototype.updateFromFirebase = function () {};

  return UnicornHat;
}(Hat_1.Hat);

exports.UnicornHat = UnicornHat;
},{"./Hat":"components/Hat.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/unicorn-hat.png":"img/unicorn-hat.png"}],"img/crickit-hat.png":[function(require,module,exports) {
module.exports = "/crickit-hat.0ec2493d.png";
},{}],"components/CrickitHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Crickit HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // @ts-ignore


var crickit_hat_png_1 = __importDefault(require("../img/crickit-hat.png"));

var CrickitHat =
/** @class */
function (_super) {
  __extends(CrickitHat, _super);

  function CrickitHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      Main_1.contextMenus.crickitHat.hat = _this;
      var menu = document.getElementById(Main_1.contextMenus.crickitHat.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      var attachMenuItem = document.getElementById("crickit-hat-attach-menu-item");
      var detachMenuItem = document.getElementById("crickit-hat-detach-menu-item");

      if (_this.raspberryPi != null) {
        attachMenuItem.className = "menu-item disabled";
        detachMenuItem.className = "menu-item";
      } else {
        var i = _this.whichRaspberryPi();

        attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
        detachMenuItem.className = "menu-item disabled";
      }
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = crickit_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  CrickitHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  CrickitHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  CrickitHat.prototype.updateFirebase = function (value) {}; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  CrickitHat.prototype.updateFromFirebase = function () {};

  return CrickitHat;
}(Hat_1.Hat);

exports.CrickitHat = CrickitHat;
},{"./Hat":"components/Hat.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/crickit-hat.png":"img/crickit-hat.png"}],"img/pan-tilt-hat.png":[function(require,module,exports) {
module.exports = "/pan-tilt-hat.76b1fdf6.png";
},{}],"components/PanTiltHat.ts":[function(require,module,exports) {
"use strict";
/*
 * Digital twin for the Unicorn HAT
 *
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Hat_1 = require("./Hat");

var Rectangle_1 = require("../math/Rectangle");

var Main_1 = require("../Main"); // @ts-ignore


var pan_tilt_hat_png_1 = __importDefault(require("../img/pan-tilt-hat.png"));

var PanTiltHat =
/** @class */
function (_super) {
  __extends(PanTiltHat, _super);

  function PanTiltHat(canvasId, uid) {
    var _this = _super.call(this, canvasId) || this;

    _this.openContextMenu = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;
      Main_1.contextMenus.panTiltHat.hat = _this;
      var menu = document.getElementById(Main_1.contextMenus.panTiltHat.id);
      menu.style.left = e.clientX + "px";
      menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
      menu.classList.add("show-menu");
      var attachMenuItem = document.getElementById("pan-tilt-hat-attach-menu-item");
      var detachMenuItem = document.getElementById("pan-tilt-hat-detach-menu-item");

      if (_this.raspberryPi != null) {
        attachMenuItem.className = "menu-item disabled";
        detachMenuItem.className = "menu-item";
      } else {
        var i = _this.whichRaspberryPi();

        attachMenuItem.className = i >= 0 ? "menu-item" : "menu-item disabled";
        detachMenuItem.className = "menu-item disabled";
      }
    };

    _this.mouseDown = function (e) {};

    _this.mouseUp = function (e) {};

    _this.mouseMove = function (e) {
      e.preventDefault();

      var rect = _this.canvas.getBoundingClientRect();

      var dx = e.clientX - rect.x;
      var dy = e.clientY - rect.y;

      if (_this.handles[0].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[0];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[1].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[1];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[2].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[2];
        _this.canvas.style.cursor = "move";
      } else if (_this.handles[3].contains(dx, dy)) {
        _this.mouseOverObject = _this.handles[3];
        _this.canvas.style.cursor = "move";
      } else {
        _this.mouseOverObject = null;
        _this.canvas.style.cursor = "default";
      }

      _this.draw();
    };

    _this.uid = uid;

    _this.canvas.addEventListener("mousedown", _this.mouseDown, false);

    _this.canvas.addEventListener("mouseup", _this.mouseUp, false);

    _this.canvas.addEventListener("mousemove", _this.mouseMove, false);

    _this.canvas.addEventListener('contextmenu', _this.openContextMenu, false);

    _this.handles.push(new Rectangle_1.Rectangle(5, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 5, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(290, 250, 30, 30));

    _this.handles.push(new Rectangle_1.Rectangle(5, 250, 30, 30));

    _this.boardImage = new Image();
    _this.boardImage.src = pan_tilt_hat_png_1.default;

    _this.setY(20);

    _this.updateFromFirebase();

    return _this;
  }

  PanTiltHat.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);
    ctx.save();
    ctx.shadowOffsetX = 8;
    ctx.shadowOffsetY = 8;
    ctx.shadowColor = "rgb(96, 96, 96)";
    ctx.shadowBlur = 8;
    ctx.drawImage(this.boardImage, 0, 0);
    ctx.restore();
    this.drawToolTips();
  };

  PanTiltHat.prototype.drawToolTips = function () {
    var context = this.canvas.getContext('2d');
    var x = 0;
    var y = -25;

    if (this.mouseOverObject == this.handles[0]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[0].getXmax() + 20;
      y += this.handles[0].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-left handle', true);
    } else if (this.mouseOverObject == this.handles[1]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[1].getXmin() - 30;
      y += this.handles[1].getYmax() + 30;
      context.drawTooltip(x, y, 20, 8, 10, 'Upper-right handle', true);
    } else if (this.mouseOverObject == this.handles[2]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[2].getXmin() - 30;
      y += this.handles[2].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-right handle', true);
    } else if (this.mouseOverObject == this.handles[3]) {
      this.drawHandle(this.mouseOverObject, context);
      x += this.handles[3].getXmax() + 20;
      y += this.handles[3].getYmin() - 5;
      context.drawTooltip(x, y, 20, 8, 10, 'Lower-left handle', true);
    }
  };

  PanTiltHat.prototype.updateFirebase = function (value) {}; // by default, sensors transmit data every second. This can be adjusted through Firebase.


  PanTiltHat.prototype.updateFromFirebase = function () {};

  return PanTiltHat;
}(Hat_1.Hat);

exports.PanTiltHat = PanTiltHat;
},{"./Hat":"components/Hat.ts","../math/Rectangle":"math/Rectangle.ts","../Main":"Main.ts","../img/pan-tilt-hat.png":"img/pan-tilt-hat.png"}],"System.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Workbench_1 = require("./Workbench");

var RaspberryPi_1 = require("./components/RaspberryPi");

var Hat_1 = require("./components/Hat");

var RainbowHat_1 = require("./components/RainbowHat");

var SenseHat_1 = require("./components/SenseHat");

var CapacitiveTouchHat_1 = require("./components/CapacitiveTouchHat");

var UnicornHat_1 = require("./components/UnicornHat");

var CrickitHat_1 = require("./components/CrickitHat");

var PanTiltHat_1 = require("./components/PanTiltHat");

var Main_1 = require("./Main");

var Rectangle_1 = require("./math/Rectangle");

var System =
/** @class */
function () {
  function System() {
    var _this = this;

    this.mcus = [];
    this.hats = [];

    this.mouseDown = function (e) {
      var rect = _this.playground.getBoundingClientRect();

      var x = e.clientX - rect.x;
      var y = e.clientY - rect.y;
      _this.selectedMovable = null; // if (this.temperatureGraph.isVisible() && this.temperatureGraph.onHandle(x - this.temperatureGraph.getX(), y - this.temperatureGraph.getY())) {
      //   this.selectedMovable = this.temperatureGraph;
      // } else if (this.pressureGraph.isVisible() && this.pressureGraph.onHandle(x - this.pressureGraph.getX(), y - this.pressureGraph.getY())) {
      //   this.selectedMovable = this.pressureGraph;
      // } else {
      // always prioritize HATs over Raspberry Pi

      for (var i = 0; i < _this.hats.length; i++) {
        if (_this.hats[i].whichHandle(x - _this.hats[i].getX(), y - _this.hats[i].getY()) >= 0) {
          _this.selectedMovable = _this.hats[i];
          break;
        }
      }

      if (_this.selectedMovable == null) {
        for (var i = 0; i < _this.mcus.length; i++) {
          if (_this.mcus[i].whichHandle(x - _this.mcus[i].getX(), y - _this.mcus[i].getY()) >= 0) {
            _this.selectedMovable = _this.mcus[i];
            break;
          }
        }
      } //}


      if (_this.selectedMovable != null) {
        _this.mouseDownRelativeX = e.clientX - _this.selectedMovable.getX();
        _this.mouseDownRelativeY = e.clientY - _this.selectedMovable.getY();
      }
    };

    this.mouseUp = function (e) {
      _this.selectedMovable = null; // close all menus upon mouse left click

      Object.keys(Main_1.contextMenus).forEach(function (key) {
        var menu = document.getElementById(Main_1.contextMenus[key].id);
        menu.classList.remove("show-menu");
      });
    };

    this.mouseLeave = function (e) {
      _this.selectedMovable = null;
    };

    this.mouseMove = function (e) {
      if (_this.selectedMovable != null) {
        _this.moveTo(e.clientX, e.clientY, _this.selectedMovable);

        _this.storeLocation(_this.selectedMovable);
      }
    };

    if (!System.database) {
      var config = {
        apiKey: "${process.env.FIREBASE_API_KEY}",
        authDomain: "${process.env.AUTH_DOMAIN}",
        projectId: "${process.env.PROJECT_ID}",
        storageBucket: "${process.env.STORAGE_BUCKET}",
        messagingSenderId: "${process.env.MESSAGING_SENDER_ID}",
        databaseURL: "https://raspberry-pi-java.firebaseio.com"
      };
      firebase.initializeApp(config); // Get a reference to the database service

      System.database = firebase.database();
    }

    this.workbench = new Workbench_1.Workbench("workbench"); //this.temperatureGraph = new LineChart("temperature-linechart", rainbowHat.temperatureSensor);
    //this.pressureGraph = new LineChart("pressure-linechart", rainbowHat.barometricPressureSensor);

    this.playground = document.getElementById("digital-twins-playground");
    this.playground.addEventListener("mousedown", this.mouseDown, false);
    this.playground.addEventListener("mouseup", this.mouseUp, false);
    this.playground.addEventListener("mousemove", this.mouseMove, false);
    document.addEventListener("mouseleave", this.mouseLeave, false); // drag and drop support

    var that = this;
    this.playground.addEventListener("dragstart", function (e) {
      that.draggedElementId = e.target.id;
    });
    this.playground.addEventListener("drag", function (e) {// console.log("dragging: " + (<HTMLElement>e.target).id);
    });
    this.playground.addEventListener("dragend", function (e) {// console.log("end drag: " + (<HTMLElement>e.target).id);
    }); // prevent default to allow drop

    this.playground.addEventListener("dragover", function (e) {
      e.preventDefault();
    }, false);
    this.playground.addEventListener("drop", function (e) {
      e.preventDefault();
      console.log("drop: " + that.draggedElementId + ", " + e.target.id);

      if (e.target.id == "workbench") {
        switch (that.draggedElementId) {
          case "raspberry-pi-image":
            that.storeLocation(that.addRaspberryPi(e.offsetX, e.offsetY, "Raspberry Pi " + Date.now().toString(16)));
            that.storeMcuSequence();
            break;

          case "rainbow-hat-image":
            that.addHatByAction("Rainbow HAT", e.offsetX, e.offsetY);
            break;

          case "sense-hat-image":
            that.addHatByAction("Sense HAT", e.offsetX, e.offsetY);
            break;

          case "capacitive-touch-hat-image":
            that.addHatByAction("Capacitive Touch HAT", e.offsetX, e.offsetY);
            break;

          case "unicorn-hat-image":
            that.addHatByAction("Unicorn HAT", e.offsetX, e.offsetY);
            break;

          case "crickit-hat-image":
            that.addHatByAction("Crickit HAT", e.offsetX, e.offsetY);
            break;

          case "pan-tilt-hat-image":
            that.addHatByAction("Pan-Tilt HAT", e.offsetX, e.offsetY);
            break;
        }
      }
    }, false);
  }

  System.prototype.addHatByAction = function (name, x, y) {
    this.storeLocation(this.addHat(name, x, y, name + " #" + Date.now().toString(16)));
    this.storeHatSequence();
  };
  /* Raspberry Pi methods */


  System.prototype.getRaspberryPiById = function (uid) {
    for (var i = 0; i < this.mcus.length; i++) {
      if (this.mcus[i] instanceof RaspberryPi_1.RaspberryPi) {
        if (this.mcus[i].uid == uid) {
          return this.mcus[i];
        }
      }
    }

    return null;
  };

  System.prototype.removeRaspberryPiByIndex = function (selectedIndex) {
    var canvas = this.mcus[selectedIndex].canvas;
    this.playground.removeChild(canvas);
    this.mcus.splice(selectedIndex, 1);
    this.storeMcuSequence();
  };

  System.prototype.removeRaspberryPi = function (raspberryPi) {
    this.removeRaspberryPiByIndex(this.mcus.indexOf(raspberryPi));
  };

  System.prototype.addRaspberryPi = function (x, y, uid) {
    var canvas = document.createElement("canvas");
    canvas.id = "raspberry-pi-" + this.mcus.length;
    canvas.width = 435;
    canvas.height = 295;
    canvas.style.display = "block";
    canvas.style.margin = "auto";
    canvas.style.position = "absolute";
    canvas.style.left = "10px";
    canvas.style.top = "10px";
    canvas.style.zIndex = "49";
    this.playground.appendChild(canvas);
    var pi = new RaspberryPi_1.RaspberryPi(canvas.id, uid);
    this.mcus.push(pi);
    pi.setX(x - canvas.width / 2);
    pi.setY(y - canvas.height / 2);
    this.draw();
    return pi;
  };

  System.prototype.whichRaspberryPi = function (x, y) {
    for (var i = 0; i < this.mcus.length; i++) {
      var mcu = this.mcus[i];

      if (mcu instanceof RaspberryPi_1.RaspberryPi) {
        var r = new Rectangle_1.Rectangle(mcu.getX(), mcu.getY(), mcu.getWidth(), mcu.getHeight());

        if (r.contains(x, y)) {
          return i;
        }
      }
    }

    return -1;
  };
  /* HAT methods */


  System.prototype.getHatById = function (uid) {
    for (var i = 0; i < this.hats.length; i++) {
      if (this.mcus[i].uid == uid) {
        return this.hats[i];
      }
    }

    return null;
  };

  System.prototype.removeHatByIndex = function (selectedIndex) {
    var canvas = this.hats[selectedIndex].canvas;
    this.playground.removeChild(canvas);
    this.hats.splice(selectedIndex, 1);
    this.storeHatSequence();
  };

  System.prototype.removeHat = function (hat) {
    this.removeHatByIndex(this.hats.indexOf(hat));
  };

  System.prototype.addHat = function (type, x, y, uid) {
    var canvas = document.createElement("canvas");
    canvas.style.display = "block";
    canvas.style.margin = "auto";
    canvas.style.position = "absolute";
    canvas.style.left = "10px";
    canvas.style.top = "10px";
    canvas.style.zIndex = "99";
    this.playground.appendChild(canvas);
    var hat = null;

    switch (type) {
      case "Rainbow HAT":
        canvas.id = "rainbow-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 290;
        hat = new RainbowHat_1.RainbowHat(canvas.id, uid);
        break;

      case "Sense HAT":
        canvas.id = "sense-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 290;
        hat = new SenseHat_1.SenseHat(canvas.id, uid);
        break;

      case "Capacitive Touch HAT":
        canvas.id = "capacitive-touch-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 290;
        hat = new CapacitiveTouchHat_1.CapacitiveTouchHat(canvas.id, uid);
        break;

      case "Unicorn HAT":
        canvas.id = "unicorn-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 380;
        hat = new UnicornHat_1.UnicornHat(canvas.id, uid);
        break;

      case "Crickit HAT":
        canvas.id = "crickit-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 290;
        hat = new CrickitHat_1.CrickitHat(canvas.id, uid);
        break;

      case "Pan-Tilt HAT":
        canvas.id = "pan-tilt-hat-" + this.hats.length;
        canvas.width = 330;
        canvas.height = 290;
        hat = new PanTiltHat_1.PanTiltHat(canvas.id, uid);
        break;
    }

    if (hat != null) {
      this.hats.push(hat);
      hat.setX(x - canvas.width / 2);
      hat.setY(y - canvas.height / 2);
      this.draw();
    }

    return hat;
  };

  System.prototype.whichHat = function (x, y) {
    for (var i = 0; i < this.hats.length; i++) {
      var hat = this.hats[i];
      var r = new Rectangle_1.Rectangle(hat.getX(), hat.getY(), hat.getWidth(), hat.getHeight());

      if (r.contains(x, y)) {
        return i;
      }
    }

    return -1;
  };

  System.prototype.draw = function () {
    this.workbench.draw();
    var i;

    for (i = 0; i < this.mcus.length; i++) {
      this.mcus[i].draw();
    }

    for (i = 0; i < this.hats.length; i++) {
      this.hats[i].draw();
    }
  };

  System.prototype.storeMcuSequence = function () {
    var s = "";

    for (var i = 0; i < this.mcus.length; i++) {
      s += this.mcus[i].getUid() + ", ";
    }

    localStorage.setItem("MCU Sequence", s.substring(0, s.length - 2));
  };

  System.prototype.storeHatSequence = function () {
    var s = "";

    for (var i = 0; i < this.hats.length; i++) {
      s += this.hats[i].getUid() + ", ";
    }

    localStorage.setItem("HAT Sequence", s.substring(0, s.length - 2));
  };

  System.prototype.storeLocation = function (m) {
    localStorage.setItem("X: " + m.getUid(), m.getX().toString());
    localStorage.setItem("Y: " + m.getUid(), m.getY().toString());

    if (m instanceof Hat_1.Hat) {
      if (m.raspberryPi != null) {
        localStorage.setItem("X: " + m.raspberryPi.getUid(), m.raspberryPi.getX().toString());
        localStorage.setItem("Y: " + m.raspberryPi.getUid(), m.raspberryPi.getY().toString());
      }
    } else if (m instanceof RaspberryPi_1.RaspberryPi) {
      if (m.hat != null) {
        localStorage.setItem("X: " + m.hat.getUid(), m.hat.getX().toString());
        localStorage.setItem("Y: " + m.hat.getUid(), m.hat.getY().toString());
      }
    }
  };

  System.prototype.moveTo = function (x, y, m) {
    var dx = x - this.mouseDownRelativeX;
    var dy = y - this.mouseDownRelativeY;
    var xmax = this.workbench.getX() + this.workbench.getWidth() - m.getWidth();

    if (dx < this.workbench.getX()) {
      dx = this.workbench.getX();
    } else if (dx > xmax) {
      dx = xmax;
    }

    var ymax = this.workbench.getY() + this.workbench.getHeight() - m.getHeight();

    if (dy < this.workbench.getY()) {
      dy = this.workbench.getY();
    } else if (dy > ymax) {
      dy = ymax;
    }

    m.setX(dx);
    m.setY(dy);

    if (m instanceof Hat_1.Hat) {
      if (m.raspberryPi != null) {
        m.raspberryPi.setX(m.getX());
        m.raspberryPi.setY(m.getY());
      }
    } else if (m instanceof RaspberryPi_1.RaspberryPi) {
      if (m.hat != null) {
        m.hat.setX(m.getX());
        m.hat.setY(m.getY());
      }
    }
  };

  return System;
}();

exports.System = System;
},{"./Workbench":"Workbench.ts","./components/RaspberryPi":"components/RaspberryPi.ts","./components/Hat":"components/Hat.ts","./components/RainbowHat":"components/RainbowHat.ts","./components/SenseHat":"components/SenseHat.ts","./components/CapacitiveTouchHat":"components/CapacitiveTouchHat.ts","./components/UnicornHat":"components/UnicornHat.ts","./components/CrickitHat":"components/CrickitHat.ts","./components/PanTiltHat":"components/PanTiltHat.ts","./Main":"Main.ts","./math/Rectangle":"math/Rectangle.ts"}],"img/full-breadboard.png":[function(require,module,exports) {
module.exports = "/full-breadboard.fbf8f80b.png";
},{}],"img/half-breadboard.png":[function(require,module,exports) {
module.exports = "/half-breadboard.d3670523.png";
},{}],"img/red-led-light.png":[function(require,module,exports) {
module.exports = "/red-led-light.57f386e5.png";
},{}],"img/green-led-light.png":[function(require,module,exports) {
module.exports = "/green-led-light.5653b4d8.png";
},{}],"img/blue-led-light.png":[function(require,module,exports) {
module.exports = "/blue-led-light.d76a5e1a.png";
},{}],"img/tricolor-led-light.png":[function(require,module,exports) {
module.exports = "/tricolor-led-light.9d522bcd.png";
},{}],"img/momentary-button.png":[function(require,module,exports) {
module.exports = "/momentary-button.d6929912.png";
},{}],"img/toggle-switch.png":[function(require,module,exports) {
module.exports = "/toggle-switch.1e723502.png";
},{}],"img/piezo-buzzer.png":[function(require,module,exports) {
module.exports = "/piezo-buzzer.b4e718d6.png";
},{}],"ComponentsPanel.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __importDefault = this && this.__importDefault || function (mod) {
  return mod && mod.__esModule ? mod : {
    "default": mod
  };
};

Object.defineProperty(exports, "__esModule", {
  value: true
}); // @ts-ignore

var raspberry_pi_png_1 = __importDefault(require("./img/raspberry-pi.png")); // @ts-ignore


var rainbow_hat_png_1 = __importDefault(require("./img/rainbow-hat.png")); // @ts-ignore


var sense_hat_png_1 = __importDefault(require("./img/sense-hat.png")); // @ts-ignore


var unicorn_hat_png_1 = __importDefault(require("./img/unicorn-hat.png")); // @ts-ignore


var crickit_hat_png_1 = __importDefault(require("./img/crickit-hat.png")); // @ts-ignore


var capacitive_touch_hat_png_1 = __importDefault(require("./img/capacitive-touch-hat.png")); // @ts-ignore


var pan_tilt_hat_png_1 = __importDefault(require("./img/pan-tilt-hat.png")); // @ts-ignore


var full_breadboard_png_1 = __importDefault(require("./img/full-breadboard.png")); // @ts-ignore


var half_breadboard_png_1 = __importDefault(require("./img/half-breadboard.png")); // @ts-ignore


var red_led_light_png_1 = __importDefault(require("./img/red-led-light.png")); // @ts-ignore


var green_led_light_png_1 = __importDefault(require("./img/green-led-light.png")); // @ts-ignore


var blue_led_light_png_1 = __importDefault(require("./img/blue-led-light.png")); // @ts-ignore


var tricolor_led_light_png_1 = __importDefault(require("./img/tricolor-led-light.png")); // @ts-ignore


var momentary_button_png_1 = __importDefault(require("./img/momentary-button.png")); // @ts-ignore


var toggle_switch_png_1 = __importDefault(require("./img/toggle-switch.png")); // @ts-ignore


var piezo_buzzer_png_1 = __importDefault(require("./img/piezo-buzzer.png"));

var ComponentsPanel =
/** @class */
function () {
  function ComponentsPanel() {}

  ComponentsPanel.prototype.getUi = function () {
    return "<h2 style=\"text-align: center; vertical-align: top; margin-top: 0;\">\n                <span style=\"font-size: 1.2em; color: teal; vertical-align: middle;\"><i class=\"fas fa-cubes\"></i></span> Components</h2>\n            <hr>\n            <div id=\"components-scroller\" style=\"overflow-y: auto; height: 360px;\">\n              <h3 style=\"text-align: left\"><span style=\"font-size: 1.2em; color: teal; vertical-align: middle\"><i class=\"fas fa-cube\"></i></span> Microcontrollers</h3>\n              <div class=\"row\" style=\"margin-right: 10px; background-color: lightskyblue; border: 1px solid #b81900; border-radius: 4px\">\n                <div class=\"column\">\n                  <img src=\"" + raspberry_pi_png_1.default + "\" draggable=\"true\" id=\"raspberry-pi-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Raspberry Pi\" alt=\"Raspberry Pi\">\n                </div>\n              </div>\n              <div class=\"vertical-divider\"></div>\n              <h3 style=\"text-align: left\"><span style=\"font-size: 1.2em; color: teal; vertical-align: middle\"><i class=\"fas fa-cube\"></i></span> HATs</h3>\n              <div class=\"row\" style=\"margin-right: 10px;  background-color: lightblue; border: 1px solid #b81900; border-radius: 4px\">\n                <div class=\"column\">\n                  <img src=\"" + rainbow_hat_png_1.default + "\" draggable=\"true\" id=\"rainbow-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Rainbow HAT\" alt=\"Rainbow HAT\">\n                  <img src=\"" + sense_hat_png_1.default + "\" draggable=\"true\" id=\"sense-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Sense HAT\" alt=\"Sense HAT\">\n                  <img src=\"" + capacitive_touch_hat_png_1.default + "\" draggable=\"true\" id=\"capacitive-touch-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Capacitive Touch HAT\" alt=\"Capacitive Touch HAT\">\n               </div>\n                <div class=\"column\">\n                  <img src=\"" + unicorn_hat_png_1.default + "\" draggable=\"true\" id=\"unicorn-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Unicorn HAT\" alt=\"Unicorn HAT\">\n                  <img src=\"" + crickit_hat_png_1.default + "\" draggable=\"true\" id=\"crickit-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Crickit HAT\" alt=\"Crickit HAT\">\n                  <img src=\"" + pan_tilt_hat_png_1.default + "\" draggable=\"true\" id=\"pan-tilt-hat-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Pan-Tilt HAT\" alt=\"Pan-Tilt HAT\">\n               </div>\n              </div>\n              <div class=\"vertical-divider\"></div>\n              <h3 style=\"text-align: left\"><span style=\"font-size: 1.2em; color: teal; vertical-align: middle\"><i class=\"fas fa-cube\"></i></span> Others</h3>\n              <div class=\"row\" style=\"margin-right: 10px;  background-color: lightyellow; border: 1px solid #b81900; border-radius: 4px\">\n                <div class=\"column\">\n                  <img src=\"" + full_breadboard_png_1.default + "\" draggable=\"true\" id=\"full-breadboard-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Full Breadboard\" alt=\"Full Breadboard\">\n                  <img src=\"" + red_led_light_png_1.default + "\" draggable=\"true\" id=\"red-led-light-image\" style=\"width:100%; cursor: pointer;\" title=\"Red LED Light\" alt=\"Red LED Light\">\n                  <img src=\"" + green_led_light_png_1.default + "\" draggable=\"true\" id=\"green-led-light-image\" style=\"width:100%; cursor: pointer;\" title=\"Green LED Light\" alt=\"Green LED Light\">\n                  <img src=\"" + blue_led_light_png_1.default + "\" draggable=\"true\" id=\"blue-led-light-image\" style=\"width:100%; cursor: pointer;\" title=\"Blue LED Light\" alt=\"Blue LED Light\">\n                 <img src=\"" + tricolor_led_light_png_1.default + "\" draggable=\"true\" id=\"tricolor-led-light-image\" style=\"width:100%; cursor: pointer;\" title=\"Tricolor LED Light\" alt=\"Tricolor LED Light\">\n                </div>\n                <div class=\"column\">\n                  <img src=\"" + half_breadboard_png_1.default + "\" draggable=\"true\" id=\"half-breadboard-image\" style=\"width:100%; cursor: pointer; box-shadow: 5px 5px 5px gray;\" title=\"Half Breadboard\" alt=\"Half Breadboard\">\n                  <img src=\"" + momentary_button_png_1.default + "\" draggable=\"true\" id=\"momentary-button-image\" style=\"width:100%; cursor: pointer;\" title=\"Momentary Button\" alt=\"Momentary Button\">\n                  <img src=\"" + toggle_switch_png_1.default + "\" draggable=\"true\" id=\"toggle-switch-image\" style=\"width:100%; cursor: pointer;\" title=\"Toggle Switch\" alt=\"Toggle Switch\">\n                  <img src=\"" + piezo_buzzer_png_1.default + "\" draggable=\"true\" id=\"piezo-buzzer-image\" style=\"width:100%; cursor: pointer;\" title=\"Piezo Buzzer\" alt=\"Piezo Buzzer\">\n                </div>\n              </div>\n            </div>";
  };

  ComponentsPanel.prototype.render = function (selectorId) {
    var element = document.getElementById(selectorId);
    element.innerHTML = this.getUi();
  };

  return ComponentsPanel;
}();

exports.ComponentsPanel = ComponentsPanel;
},{"./img/raspberry-pi.png":"img/raspberry-pi.png","./img/rainbow-hat.png":"img/rainbow-hat.png","./img/sense-hat.png":"img/sense-hat.png","./img/unicorn-hat.png":"img/unicorn-hat.png","./img/crickit-hat.png":"img/crickit-hat.png","./img/capacitive-touch-hat.png":"img/capacitive-touch-hat.png","./img/pan-tilt-hat.png":"img/pan-tilt-hat.png","./img/full-breadboard.png":"img/full-breadboard.png","./img/half-breadboard.png":"img/half-breadboard.png","./img/red-led-light.png":"img/red-led-light.png","./img/green-led-light.png":"img/green-led-light.png","./img/blue-led-light.png":"img/blue-led-light.png","./img/tricolor-led-light.png":"img/tricolor-led-light.png","./img/momentary-button.png":"img/momentary-button.png","./img/toggle-switch.png":"img/toggle-switch.png","./img/piezo-buzzer.png":"img/piezo-buzzer.png"}],"ui/ComponentContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var ComponentContextMenu =
/** @class */
function () {
  function ComponentContextMenu() {}

  ComponentContextMenu.prototype.render = function (selectorId) {
    var element = document.getElementById(selectorId);
    element.innerHTML = this.getUi();
  };

  return ComponentContextMenu;
}();

exports.ComponentContextMenu = ComponentContextMenu;
},{}],"ui/HatContextMenu.ts":[function(require,module,exports) {
"use strict";

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Main_1 = require("../Main");

var ComponentContextMenu_1 = require("./ComponentContextMenu");

var HatContextMenu =
/** @class */
function (_super) {
  __extends(HatContextMenu, _super);

  function HatContextMenu() {
    return _super !== null && _super.apply(this, arguments) || this;
  }

  HatContextMenu.prototype.addListeners = function () {
    var attachButton = document.getElementById(this.id + "-attach-button");
    attachButton.addEventListener("click", this.attachButtonClick.bind(this), false);
    var detachButton = document.getElementById(this.id + "-detach-button");
    detachButton.addEventListener("click", this.detachButtonClick.bind(this), false);
    var deleteButton = document.getElementById(this.id + "-delete-button");
    deleteButton.addEventListener("click", this.deleteButtonClick.bind(this), false);
  };

  HatContextMenu.prototype.deleteButtonClick = function (e) {
    if (this.hat) {
      if (confirm("Are you sure you want to delete " + this.hat.uid + "?")) {
        Main_1.system.removeHat(this.hat);
      }
    }
  };

  HatContextMenu.prototype.attachButtonClick = function (e) {
    e.preventDefault();
    var menu = document.getElementById(this.id);
    menu.classList.remove("show-menu");
    var i = this.hat.whichRaspberryPi();

    if (i >= 0) {
      this.hat.attach(Main_1.system.mcus[i]);
    }
  };

  HatContextMenu.prototype.detachButtonClick = function (e) {
    e.preventDefault();
    var menu = document.getElementById(this.id);
    menu.classList.remove("show-menu");
    this.hat.attach(null);
  };

  return HatContextMenu;
}(ComponentContextMenu_1.ComponentContextMenu);

exports.HatContextMenu = HatContextMenu;
},{"../Main":"Main.ts","./ComponentContextMenu":"ui/ComponentContextMenu.ts"}],"ui/RainbowHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var RainbowHatContextMenu =
/** @class */
function (_super) {
  __extends(RainbowHatContextMenu, _super);

  function RainbowHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "rainbow-hat-context-menu";
    return _this;
  }

  RainbowHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"rainbow-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"rainbow-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"rainbow-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return RainbowHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.RainbowHatContextMenu = RainbowHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"ui/WorkbenchContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var ComponentContextMenu_1 = require("./ComponentContextMenu");

var WorkbenchContextMenu =
/** @class */
function (_super) {
  __extends(WorkbenchContextMenu, _super);

  function WorkbenchContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "workbench-context-menu";
    return _this;
  }

  WorkbenchContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-folder-open\"></i><span class=\"menu-text\">Open</span></button>\n              </li>\n              <li class=\"menu-item disabled\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-download\"></i><span class=\"menu-text\">Save</span></button>\n              </li>\n              <li class=\"menu-separator\"></li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item submenu\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-file-import\"></i><span class=\"menu-text\">Import</span></button>\n\n                <menu class=\"menu\" style=\"width: 160px;\">\n                  <li class=\"menu-item\">\n                    <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Breadboard</span></button>\n                  </li>\n\n                  <li class=\"menu-item submenu\">\n                    <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Sensors</span></button>\n                    <menu class=\"menu\" style=\"width: 120px;\">\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">BME280</span></button>\n                      </li>\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">TSL2561</span></button>\n                      </li>\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">HCSR04</span></button>\n                      </li>\n                    </menu>\n                  </li>\n\n                  <li class=\"menu-item submenu\">\n                    <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Actuators</span></button>\n                    <menu class=\"menu\" style=\"width: 180px;\">\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Buzzer</span></button>\n                      </li>\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Servo Motor</span></button>\n                      </li>\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">LED Light</span></button>\n                      </li>\n                      <li class=\"menu-item\">\n                        <button type=\"button\" class=\"menu-btn\"><span class=\"menu-text\">Multicolor LED Light</span></button>\n                      </li>\n                    </menu>\n                  </li>\n\n                </menu>\n              </li>\n            </menu>";
  };

  WorkbenchContextMenu.prototype.addListeners = function () {// TODO
  };

  return WorkbenchContextMenu;
}(ComponentContextMenu_1.ComponentContextMenu);

exports.WorkbenchContextMenu = WorkbenchContextMenu;
},{"./ComponentContextMenu":"ui/ComponentContextMenu.ts"}],"ui/LineChartContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var LineChartContextMenu =
/** @class */
function () {
  function LineChartContextMenu() {
    this.id = "linechart-context-menu";
  }

  LineChartContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"z-index: 10000\">\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-palette\"></i><span class=\"menu-text\">Colors</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\"><i class=\"fas fa-cog\"></i><span class=\"menu-text\">Options</span></button>\n              </li>\n            </menu>";
  };

  LineChartContextMenu.prototype.render = function (selectorId) {
    var element = document.getElementById(selectorId);
    element.innerHTML = this.getUi();
  };

  return LineChartContextMenu;
}();

exports.LineChartContextMenu = LineChartContextMenu;
},{}],"ui/RaspberryPiContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Main_1 = require("../Main");

var ComponentContextMenu_1 = require("./ComponentContextMenu");

var RaspberryPiContextMenu =
/** @class */
function (_super) {
  __extends(RaspberryPiContextMenu, _super);

  function RaspberryPiContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "raspberry-pi-context-menu";
    return _this;
  }

  RaspberryPiContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\">\n                <button type=\"button\" id=\"raspberry-pi-context-menu-delete-button\" class=\"menu-btn\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" id=\"raspberry-pi-context-menu-settings-button\" class=\"menu-btn\"><i class=\"fas fa-cog\"></i><span class=\"menu-text\">Settings</span></button>\n              </li>\n            </menu>";
  };

  RaspberryPiContextMenu.prototype.addListeners = function () {
    var deleteButton = document.getElementById("raspberry-pi-context-menu-delete-button");
    deleteButton.addEventListener("click", this.deleteButtonClick.bind(this), false);
  };

  RaspberryPiContextMenu.prototype.deleteButtonClick = function (e) {
    if (this.raspberryPi) {
      if (confirm("Are you sure you want to delete " + this.raspberryPi.uid + "?")) {
        Main_1.system.removeRaspberryPi(this.raspberryPi);
      }
    }
  };

  return RaspberryPiContextMenu;
}(ComponentContextMenu_1.ComponentContextMenu);

exports.RaspberryPiContextMenu = RaspberryPiContextMenu;
},{"../Main":"Main.ts","./ComponentContextMenu":"ui/ComponentContextMenu.ts"}],"ui/ColorPickerContextMenu.ts":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
/*
 * @author Charles Xie
 */

var Main_1 = require("../Main");

var ColorPickerContextMenu =
/** @class */
function () {
  function ColorPickerContextMenu() {
    this.id = "colorpicker-context-menu";
  }

  ColorPickerContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 338px; z-index: 10000;\">\n              <li class=\"menu-item\">\n                <div id=\"colorpicker-title\" style=\"text-align: center; padding: 4px; font-family: inherit; font-size: 14px;\"></div>\n                <div id=\"colorpicker\" style=\"cursor: crosshair; margin: 1px 1px 1px 1px\">\n                  <canvas id=\"color-block\" height=\"300\" width=\"300\"></canvas>\n                  <canvas id=\"color-strip\" height=\"300\" width=\"30\"></canvas>\n                </div>\n                <div style=\"display: table; margin: auto; padding: 5px 5px 5px 5px;\">\n                  <input type=\"text\" readonly id=\"colorpicker-hex-code\" value=\"#FFFFFF\" style=\"width: 60px; height: 20px; border: 1px solid black; vertical-align: middle; text-align: center; font: 12px Verdana;\">\n                  <div class=\"horizontal-divider\"></div>\n                  <div id=\"colorpicker-label\" style=\"vertical-align: middle; width: 20px; height: 20px; border: 2px solid black; display: inline-block;\"></div>\n                  <div class=\"horizontal-divider\"></div>\n                  <button id=\"colorpicker-cancel-button\" style=\"font: 12px Verdana\">Cancel</button>\n                  <div class=\"horizontal-divider\"></div>\n                  <button id=\"colorpicker-ok-button\" style=\"font: 12px Verdana\">OK</button>\n                </div>\n              </li>\n            </menu>";
  };

  ColorPickerContextMenu.prototype.render = function (selectorId) {
    var element = document.getElementById(selectorId);
    element.innerHTML = this.getUi();
    var cancelButton = document.getElementById("colorpicker-cancel-button");

    cancelButton.onclick = function () {
      var menu = document.getElementById(that.id);
      menu.classList.remove("show-menu");
    };

    var okButton = document.getElementById("colorpicker-ok-button");
    var that = this;

    okButton.onclick = function () {
      if (that.rainbowHat) {
        that.rainbowHat.setSelectedRgbLedLightColor(Main_1.system.colorPicker.getSelectedColor());
      }

      var menu = document.getElementById(that.id);
      menu.classList.remove("show-menu");
    }; // prevent the mouse event from being propagated to the playground


    var menu = document.getElementById(this.id);

    menu.onmousedown = function (e) {
      e.stopPropagation();
    };

    menu.onmouseup = function (e) {
      e.stopPropagation();
    };

    menu.onmousemove = function (e) {
      e.stopPropagation();
    };
  };

  return ColorPickerContextMenu;
}();

exports.ColorPickerContextMenu = ColorPickerContextMenu;
},{"../Main":"Main.ts"}],"code/Block.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Block =
/** @class */
function () {
  function Block(x, y, width, height, name) {
    this.radius = 10;
    this.margin = 30;
    this.x = x;
    this.y = y;
    this.width = width;
    this.height = height;
    this.name = name;
  }

  Block.prototype.draw = function (ctx) {
    ctx.fillStyle = "lightgray";
    ctx.fillRoundedRect(this.x, this.y, this.width, this.height, this.radius);
    ctx.lineWidth = 1;
    ctx.strokeStyle = "black";
    ctx.drawRoundedRect(this.x, this.y, this.width, this.height, this.radius);
    ctx.fillStyle = "white";
    ctx.beginPath();
    ctx.rect(this.x + this.margin, this.y + this.margin, this.width - 2 * this.margin, this.height - 2 * this.margin);
    ctx.closePath();
    ctx.fill();
    ctx.strokeStyle = "black";
    ctx.stroke();
    ctx.fillStyle = "black";
    ctx.font = "20px Arial";
    var textMetrics = ctx.measureText(this.name);
    ctx.fillText(this.name, this.x + this.width / 2 - textMetrics.width / 2, this.y + this.height / 2 + 8);
    var r = 8;
    ctx.beginPath();
    ctx.arc(this.x, this.y + this.height / 3, r, 0.5 * Math.PI, 1.5 * Math.PI, false);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(this.x, this.y + this.height * 2 / 3, r, 0.5 * Math.PI, 1.5 * Math.PI, false);
    ctx.stroke();
    ctx.beginPath();
    ctx.arc(this.x + this.width, this.y + this.height / 2, r, 0.5 * Math.PI, 1.5 * Math.PI, true);
    ctx.stroke();
  };

  return Block;
}();

exports.Block = Block;
},{}],"math/Point.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Point =
/** @class */
function () {
  function Point(x, y) {
    this.x = x;
    this.y = y;
  }

  Point.prototype.toString = function () {
    return "(" + this.x + ", " + this.y + ")";
  };

  return Point;
}();

exports.Point = Point;
},{}],"code/Codespace.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

Object.defineProperty(exports, "__esModule", {
  value: true
});

var Block_1 = require("./Block");

var Point_1 = require("../math/Point");

var Codespace =
/** @class */
function () {
  function Codespace(canvasId) {
    this.blocks = [];
    this.canvas = document.getElementById(canvasId);
    this.canvas.addEventListener("mousedown", this.mouseDown.bind(this), false);
    this.canvas.addEventListener("mouseup", this.mouseUp.bind(this), false);
    this.canvas.addEventListener("mousemove", this.mouseMove.bind(this), false);
    this.canvas.addEventListener('contextmenu', this.openContextMenu.bind(this), false);
    this.blocks.push(new Block_1.Block(20, 20, 160, 100, "X + Y"));
    this.blocks.push(new Block_1.Block(220, 220, 160, 100, "X * Y"));
  }

  Codespace.prototype.draw = function () {
    var ctx = this.canvas.getContext('2d');
    ctx.clearRect(0, 0, this.canvas.width, this.canvas.height);

    for (var i = 0; i < this.blocks.length; i++) {
      this.blocks[i].draw(ctx);
    }

    var points = [];
    points.push(new Point_1.Point(188, 70));
    points.push(new Point_1.Point(240, 140));
    points.push(new Point_1.Point(20, 200));
    points.push(new Point_1.Point(212, 288));
    this.drawSpline(points, ctx);
  };

  Codespace.prototype.drawSpline = function (points, ctx) {
    ctx.beginPath();
    ctx.moveTo(points[0].x, points[0].y);
    var i;

    for (i = 1; i < points.length - 2; i++) {
      var xc = (points[i].x + points[i + 1].x) / 2;
      var yc = (points[i].y + points[i + 1].y) / 2;
      ctx.quadraticCurveTo(points[i].x, points[i].y, xc, yc);
    }

    ctx.quadraticCurveTo(points[i].x, points[i].y, points[i + 1].x, points[i + 1].y);
    ctx.stroke();
  }; // detect if (x, y) is inside this workbench


  Codespace.prototype.contains = function (x, y) {
    return x > this.canvas.offsetLeft && x < this.canvas.offsetLeft + this.canvas.width && y > this.canvas.offsetTop && y < this.canvas.offsetTop + this.canvas.height;
  };

  Codespace.prototype.getX = function () {
    return 10;
  };

  Codespace.prototype.getY = function () {
    return 10;
  };

  Codespace.prototype.getWidth = function () {
    return this.canvas.width;
  };

  Codespace.prototype.getHeight = function () {
    return this.canvas.height;
  };

  Codespace.prototype.mouseDown = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
  };

  Codespace.prototype.mouseUp = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
    var context = this.canvas.getContext("2d");
  };

  Codespace.prototype.mouseMove = function (e) {
    e.preventDefault();
    var rect = this.canvas.getBoundingClientRect();
    var dx = e.clientX - rect.x;
    var dy = e.clientY - rect.y;
    var context = this.canvas.getContext("2d");
  };

  Codespace.prototype.openContextMenu = function (e) {
    e.preventDefault();
    var menu = document.getElementById("workbench-context-menu");
    menu.style.left = e.clientX + "px";
    menu.style.top = e.clientY - document.getElementById("tabs").getBoundingClientRect().bottom + "px";
    menu.classList.add("show-menu");
  };

  return Codespace;
}();

exports.Codespace = Codespace;
},{"./Block":"code/Block.ts","../math/Point":"math/Point.ts"}],"code/Code.ts":[function(require,module,exports) {
"use strict";

Object.defineProperty(exports, "__esModule", {
  value: true
});
/*
 * @author Charles Xie
 */

var Codespace_1 = require("./Codespace");

var Code =
/** @class */
function () {
  function Code() {
    this.codespace = new Codespace_1.Codespace("codespace");
  }

  Code.prototype.draw = function () {
    this.codespace.draw();
  };

  return Code;
}();

exports.Code = Code;
},{"./Codespace":"code/Codespace.ts"}],"ui/SenseHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var SenseHatContextMenu =
/** @class */
function (_super) {
  __extends(SenseHatContextMenu, _super);

  function SenseHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "sense-hat-context-menu";
    return _this;
  }

  SenseHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"sense-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"sense-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"sense-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return SenseHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.SenseHatContextMenu = SenseHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"ui/CapacitiveTouchHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var CapacitiveTouchHatContextMenu =
/** @class */
function (_super) {
  __extends(CapacitiveTouchHatContextMenu, _super);

  function CapacitiveTouchHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "capacitive-touch-hat-context-menu";
    return _this;
  }

  CapacitiveTouchHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"capacitive-touch-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"capacitive-touch-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"capacitive-touch-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return CapacitiveTouchHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.CapacitiveTouchHatContextMenu = CapacitiveTouchHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"ui/UnicornHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var UnicornHatContextMenu =
/** @class */
function (_super) {
  __extends(UnicornHatContextMenu, _super);

  function UnicornHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "unicorn-hat-context-menu";
    return _this;
  }

  UnicornHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"unicorn-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"unicorn-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"unicorn-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return UnicornHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.UnicornHatContextMenu = UnicornHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"ui/CrickitHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var CrickitHatContextMenu =
/** @class */
function (_super) {
  __extends(CrickitHatContextMenu, _super);

  function CrickitHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "crickit-hat-context-menu";
    return _this;
  }

  CrickitHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"crickit-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"crickit-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"crickit-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return CrickitHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.CrickitHatContextMenu = CrickitHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"ui/PanTiltHatContextMenu.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __extends = this && this.__extends || function () {
  var _extendStatics = function extendStatics(d, b) {
    _extendStatics = Object.setPrototypeOf || {
      __proto__: []
    } instanceof Array && function (d, b) {
      d.__proto__ = b;
    } || function (d, b) {
      for (var p in b) {
        if (b.hasOwnProperty(p)) d[p] = b[p];
      }
    };

    return _extendStatics(d, b);
  };

  return function (d, b) {
    _extendStatics(d, b);

    function __() {
      this.constructor = d;
    }

    d.prototype = b === null ? Object.create(b) : (__.prototype = b.prototype, new __());
  };
}();

Object.defineProperty(exports, "__esModule", {
  value: true
});

var HatContextMenu_1 = require("./HatContextMenu");

var PanTiltHatContextMenu =
/** @class */
function (_super) {
  __extends(PanTiltHatContextMenu, _super);

  function PanTiltHatContextMenu() {
    var _this = _super.call(this) || this;

    _this.id = "pan-tilt-hat-context-menu";
    return _this;
  }

  PanTiltHatContextMenu.prototype.getUi = function () {
    return "<menu id=\"" + this.id + "\" class=\"menu\" style=\"width: 120px; z-index: 10000\">\n              <li class=\"menu-item\" id=\"pan-tilt-hat-attach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-attach-button\"><i class=\"fas fa-angle-double-down\"></i><span class=\"menu-text\">Attach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"pan-tilt-hat-detach-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-detach-button\"><i class=\"fas fa-angle-double-up\"></i><span class=\"menu-text\">Detach</span></button>\n              </li>\n              <li class=\"menu-item\" id=\"pan-tilt-hat-delete-menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-delete-button\"><i class=\"fas fa-trash\"></i><span class=\"menu-text\">Delete</span></button>\n              </li>\n              <li class=\"menu-item\">\n                <button type=\"button\" class=\"menu-btn\" id=\"" + this.id + "-code-button\"><i class=\"fas fa-code\"></i><span class=\"menu-text\">Code</span></button>\n              </li>\n            </menu>";
  };

  return PanTiltHatContextMenu;
}(HatContextMenu_1.HatContextMenu);

exports.PanTiltHatContextMenu = PanTiltHatContextMenu;
},{"./HatContextMenu":"ui/HatContextMenu.ts"}],"Main.ts":[function(require,module,exports) {
"use strict";
/*
 * @author Charles Xie
 */

var __importStar = this && this.__importStar || function (mod) {
  if (mod && mod.__esModule) return mod;
  var result = {};
  if (mod != null) for (var k in mod) {
    if (Object.hasOwnProperty.call(mod, k)) result[k] = mod[k];
  }
  result["default"] = mod;
  return result;
};

Object.defineProperty(exports, "__esModule", {
  value: true
});

require("@fortawesome/fontawesome-free/css/all.css");

var Constants = __importStar(require("./Constants"));

var User_1 = require("./User");

var System_1 = require("./System");

var ComponentsPanel_1 = require("./ComponentsPanel");

var RainbowHatContextMenu_1 = require("./ui/RainbowHatContextMenu");

var WorkbenchContextMenu_1 = require("./ui/WorkbenchContextMenu");

var LineChartContextMenu_1 = require("./ui/LineChartContextMenu");

var RaspberryPiContextMenu_1 = require("./ui/RaspberryPiContextMenu");

var ColorPickerContextMenu_1 = require("./ui/ColorPickerContextMenu");

var Code_1 = require("./code/Code");

var SenseHatContextMenu_1 = require("./ui/SenseHatContextMenu");

var CapacitiveTouchHatContextMenu_1 = require("./ui/CapacitiveTouchHatContextMenu");

var UnicornHatContextMenu_1 = require("./ui/UnicornHatContextMenu");

var CrickitHatContextMenu_1 = require("./ui/CrickitHatContextMenu");

var PanTiltHatContextMenu_1 = require("./ui/PanTiltHatContextMenu");

exports.system = new System_1.System();
exports.code = new Code_1.Code();
exports.user = new User_1.User("Charles", null, "Xie");
exports.contextMenus = {};
var social = "<span style=\"font-size: 2em; vertical-align: middle; cursor: pointer;\"><i class=\"fab fa-facebook-square\"></i></span>\n              <span style=\"font-size: 2em; vertical-align: middle; cursor: pointer;\"><i class=\"fab fa-weixin\"></i></span>\n              <span style=\"font-size: 2em; vertical-align: middle; cursor: pointer;\"><i class=\"fab fa-twitter\"></i></span>\n              <span style=\"font-size: 2em; vertical-align: middle; cursor: pointer;\"><i class=\"fab fa-weibo\"></i></span>\n              <span style=\"font-size: 2em; vertical-align: middle; cursor: pointer;\"><i class=\"fab fa-youtube\"></i></span>";

window.onload = function () {
  var signinLabel = document.getElementById("sign-in-label");
  signinLabel.innerHTML = "Hello, " + exports.user.firstName;
  var nameLabel = document.getElementById("name-label");
  nameLabel.innerHTML = Constants.Software.name;
  var versionLabel = document.getElementById("version-label");
  versionLabel.innerHTML = Constants.Software.version;
  var creditLabel = document.getElementById('credit');
  creditLabel.innerHTML = social + "<div class='horizontal-divider'></div>" + Constants.Software.name + " " + Constants.Software.version + ", " + exports.user.fullName + " , &copy; " + new Date().getFullYear();
  var digitalTwinsTabButton = document.getElementById("digital-twins-tab-button");
  digitalTwinsTabButton.addEventListener("click", function (event) {
    selectTab(digitalTwinsTabButton, "digital-twins-playground");
  });
  var diagramTabButton = document.getElementById("diagram-tab-button");
  diagramTabButton.addEventListener("click", function (event) {
    selectTab(diagramTabButton, "diagram-playground");
  });
  var codeTabButton = document.getElementById("code-tab-button");
  codeTabButton.addEventListener("click", function (event) {
    selectTab(codeTabButton, "code-playground");
  });
  var workbenchContextMenu = new WorkbenchContextMenu_1.WorkbenchContextMenu();
  workbenchContextMenu.render("workbench-context-menu-placeholder");
  exports.contextMenus.workbench = workbenchContextMenu;
  var raspberryPiContextMenu = new RaspberryPiContextMenu_1.RaspberryPiContextMenu();
  raspberryPiContextMenu.render("raspberry-pi-context-menu-placeholder");
  raspberryPiContextMenu.addListeners();
  exports.contextMenus.raspberryPi = raspberryPiContextMenu;
  var rainbowHatContextMenu = new RainbowHatContextMenu_1.RainbowHatContextMenu();
  rainbowHatContextMenu.render("rainbow-hat-context-menu-placeholder");
  rainbowHatContextMenu.addListeners();
  exports.contextMenus.rainbowHat = rainbowHatContextMenu;
  var senseHatContextMenu = new SenseHatContextMenu_1.SenseHatContextMenu();
  senseHatContextMenu.render("sense-hat-context-menu-placeholder");
  senseHatContextMenu.addListeners();
  exports.contextMenus.senseHat = senseHatContextMenu;
  var unicornHatContextMenu = new UnicornHatContextMenu_1.UnicornHatContextMenu();
  unicornHatContextMenu.render("unicorn-hat-context-menu-placeholder");
  unicornHatContextMenu.addListeners();
  exports.contextMenus.unicornHat = unicornHatContextMenu;
  var crickitHatContextMenu = new CrickitHatContextMenu_1.CrickitHatContextMenu();
  crickitHatContextMenu.render("crickit-hat-context-menu-placeholder");
  crickitHatContextMenu.addListeners();
  exports.contextMenus.crickitHat = crickitHatContextMenu;
  var panTiltHatContextMenu = new PanTiltHatContextMenu_1.PanTiltHatContextMenu();
  panTiltHatContextMenu.render("pan-tilt-hat-context-menu-placeholder");
  panTiltHatContextMenu.addListeners();
  exports.contextMenus.panTiltHat = panTiltHatContextMenu;
  var capacitiveTouchHatContextMenu = new CapacitiveTouchHatContextMenu_1.CapacitiveTouchHatContextMenu();
  capacitiveTouchHatContextMenu.render("capacitive-touch-hat-context-menu-placeholder");
  capacitiveTouchHatContextMenu.addListeners();
  exports.contextMenus.capacitiveTouchHat = capacitiveTouchHatContextMenu;
  var lineChartContextMenu = new LineChartContextMenu_1.LineChartContextMenu();
  lineChartContextMenu.render("linechart-context-menu-placeholder");
  exports.contextMenus.lineChart = lineChartContextMenu;
  var colorPickerContextMenu = new ColorPickerContextMenu_1.ColorPickerContextMenu();
  colorPickerContextMenu.render("colorpicker-context-menu-placeholder");
  exports.contextMenus.colorPicker = colorPickerContextMenu;
  var componentsPanel = new ComponentsPanel_1.ComponentsPanel();
  componentsPanel.render("digital-twins-playground-components-panel"); // read locally stored properties

  restoreMcus();
  restoreHats(); // restoreLocation(system.temperatureGraph);
  // restoreLocation(system.pressureGraph);
  // restoreVisibility(system.temperatureGraph);
  // restoreVisibility(system.pressureGraph);

  resize();
  draw();
};

function selectTab(button, tabId) {
  // Get all elements with class="tabcontent" and hide them
  var tabcontent = document.getElementsByClassName("tabcontent");

  for (var i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
  } // Get all elements with class="tablinks" and remove the class "active"


  var tablinks = document.getElementsByClassName("tablinks");

  for (var i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  } // Show the current tab, and add an "active" class to the button that opened the tab


  document.getElementById(tabId).style.display = "block";
  button.className += " active";
}

function restoreVisibility(g) {
  var x = localStorage.getItem("Visible: " + g.getUid());

  if (x != null) {
    g.setVisible("true" == x);
    g.draw();
  }
}

function restoreHats() {
  var s = localStorage.getItem("HAT Sequence");

  if (s != null) {
    var t = s.split(",");

    if (t.length > 0) {
      exports.system.hats = [];
    }

    for (var i = 0; i < t.length; i++) {
      t[i] = t[i].trim();
      var name = t[i].substring(0, t[i].indexOf("#") - 1);
      exports.system.addHat(name, 0, 0, t[i]);
    }
  }

  restoreLocations(exports.system.hats);

  for (var i = 0; i < exports.system.hats.length; i++) {
    var id = localStorage.getItem("Attachment: " + exports.system.hats[i].getUid());

    if (id != null) {
      var pi = exports.system.getRaspberryPiById(id);

      if (pi) {
        exports.system.hats[i].attach(pi);
      }
    }
  }

  setTimeout(function () {
    exports.system.draw();
  }, 100);
}

function restoreMcus() {
  var s = localStorage.getItem("MCU Sequence");

  if (s != null) {
    var t = s.split(",");

    if (t.length > 0) {
      exports.system.mcus = [];
    }

    for (var i = 0; i < t.length; i++) {
      t[i] = t[i].trim();

      if (t[i].startsWith("Raspberry Pi")) {
        exports.system.addRaspberryPi(0, 0, t[i]);
      }
    }
  }

  restoreLocations(exports.system.mcus);
  setTimeout(function () {
    exports.system.draw();
  }, 0);
}

function restoreLocations(m) {
  for (var i = 0; i < m.length; i++) {
    restoreLocation(m[i]);
  }
}

function restoreLocation(m) {
  var x = localStorage.getItem("X: " + m.getUid());

  if (x != null) {
    m.setX(parseInt(x));
  }

  var y = localStorage.getItem("Y: " + m.getUid());

  if (y != null) {
    m.setY(parseInt(y));
  }
}

window.onresize = function () {
  resize();
  draw();
};

function resize() {
  var workbenchRect = exports.system.workbench.canvas.getBoundingClientRect();
  exports.system.workbench.canvas.width = window.innerWidth - 2 * workbenchRect.left - 4;
  exports.system.workbench.canvas.height = window.innerHeight - workbenchRect.top - 50;
  var codespaceRect = exports.code.codespace.canvas.getBoundingClientRect();
  exports.code.codespace.canvas.width = window.innerWidth - 2 * workbenchRect.left - 4;
  exports.code.codespace.canvas.height = window.innerHeight - workbenchRect.top - 50;
  var componentsScroller = document.getElementById("components-scroller");
  componentsScroller.style.height = exports.system.workbench.canvas.height * 0.8 + "px";
}

function draw() {
  exports.system.draw();
  exports.code.draw();
}
},{"@fortawesome/fontawesome-free/css/all.css":"node_modules/@fortawesome/fontawesome-free/css/all.css","./Constants":"Constants.ts","./User":"User.ts","./System":"System.ts","./ComponentsPanel":"ComponentsPanel.ts","./ui/RainbowHatContextMenu":"ui/RainbowHatContextMenu.ts","./ui/WorkbenchContextMenu":"ui/WorkbenchContextMenu.ts","./ui/LineChartContextMenu":"ui/LineChartContextMenu.ts","./ui/RaspberryPiContextMenu":"ui/RaspberryPiContextMenu.ts","./ui/ColorPickerContextMenu":"ui/ColorPickerContextMenu.ts","./code/Code":"code/Code.ts","./ui/SenseHatContextMenu":"ui/SenseHatContextMenu.ts","./ui/CapacitiveTouchHatContextMenu":"ui/CapacitiveTouchHatContextMenu.ts","./ui/UnicornHatContextMenu":"ui/UnicornHatContextMenu.ts","./ui/CrickitHatContextMenu":"ui/CrickitHatContextMenu.ts","./ui/PanTiltHatContextMenu":"ui/PanTiltHatContextMenu.ts"}],"node_modules/parcel/src/builtins/hmr-runtime.js":[function(require,module,exports) {
var global = arguments[3];
var OVERLAY_ID = '__parcel__error__overlay__';
var OldModule = module.bundle.Module;

function Module(moduleName) {
  OldModule.call(this, moduleName);
  this.hot = {
    data: module.bundle.hotData,
    _acceptCallbacks: [],
    _disposeCallbacks: [],
    accept: function (fn) {
      this._acceptCallbacks.push(fn || function () {});
    },
    dispose: function (fn) {
      this._disposeCallbacks.push(fn);
    }
  };
  module.bundle.hotData = null;
}

module.bundle.Module = Module;
var checkedAssets, assetsToAccept;
var parent = module.bundle.parent;

if ((!parent || !parent.isParcelRequire) && typeof WebSocket !== 'undefined') {
  var hostname = "" || location.hostname;
  var protocol = location.protocol === 'https:' ? 'wss' : 'ws';
  var ws = new WebSocket(protocol + '://' + hostname + ':' + "53169" + '/');

  ws.onmessage = function (event) {
    checkedAssets = {};
    assetsToAccept = [];
    var data = JSON.parse(event.data);

    if (data.type === 'update') {
      var handled = false;
      data.assets.forEach(function (asset) {
        if (!asset.isNew) {
          var didAccept = hmrAcceptCheck(global.parcelRequire, asset.id);

          if (didAccept) {
            handled = true;
          }
        }
      }); // Enable HMR for CSS by default.

      handled = handled || data.assets.every(function (asset) {
        return asset.type === 'css' && asset.generated.js;
      });

      if (handled) {
        console.clear();
        data.assets.forEach(function (asset) {
          hmrApply(global.parcelRequire, asset);
        });
        assetsToAccept.forEach(function (v) {
          hmrAcceptRun(v[0], v[1]);
        });
      } else if (location.reload) {
        // `location` global exists in a web worker context but lacks `.reload()` function.
        location.reload();
      }
    }

    if (data.type === 'reload') {
      ws.close();

      ws.onclose = function () {
        location.reload();
      };
    }

    if (data.type === 'error-resolved') {
      console.log('[parcel] ✨ Error resolved');
      removeErrorOverlay();
    }

    if (data.type === 'error') {
      console.error('[parcel] 🚨  ' + data.error.message + '\n' + data.error.stack);
      removeErrorOverlay();
      var overlay = createErrorOverlay(data);
      document.body.appendChild(overlay);
    }
  };
}

function removeErrorOverlay() {
  var overlay = document.getElementById(OVERLAY_ID);

  if (overlay) {
    overlay.remove();
  }
}

function createErrorOverlay(data) {
  var overlay = document.createElement('div');
  overlay.id = OVERLAY_ID; // html encode message and stack trace

  var message = document.createElement('div');
  var stackTrace = document.createElement('pre');
  message.innerText = data.error.message;
  stackTrace.innerText = data.error.stack;
  overlay.innerHTML = '<div style="background: black; font-size: 16px; color: white; position: fixed; height: 100%; width: 100%; top: 0px; left: 0px; padding: 30px; opacity: 0.85; font-family: Menlo, Consolas, monospace; z-index: 9999;">' + '<span style="background: red; padding: 2px 4px; border-radius: 2px;">ERROR</span>' + '<span style="top: 2px; margin-left: 5px; position: relative;">🚨</span>' + '<div style="font-size: 18px; font-weight: bold; margin-top: 20px;">' + message.innerHTML + '</div>' + '<pre>' + stackTrace.innerHTML + '</pre>' + '</div>';
  return overlay;
}

function getParents(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return [];
  }

  var parents = [];
  var k, d, dep;

  for (k in modules) {
    for (d in modules[k][1]) {
      dep = modules[k][1][d];

      if (dep === id || Array.isArray(dep) && dep[dep.length - 1] === id) {
        parents.push(k);
      }
    }
  }

  if (bundle.parent) {
    parents = parents.concat(getParents(bundle.parent, id));
  }

  return parents;
}

function hmrApply(bundle, asset) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (modules[asset.id] || !bundle.parent) {
    var fn = new Function('require', 'module', 'exports', asset.generated.js);
    asset.isNew = !modules[asset.id];
    modules[asset.id] = [fn, asset.deps];
  } else if (bundle.parent) {
    hmrApply(bundle.parent, asset);
  }
}

function hmrAcceptCheck(bundle, id) {
  var modules = bundle.modules;

  if (!modules) {
    return;
  }

  if (!modules[id] && bundle.parent) {
    return hmrAcceptCheck(bundle.parent, id);
  }

  if (checkedAssets[id]) {
    return;
  }

  checkedAssets[id] = true;
  var cached = bundle.cache[id];
  assetsToAccept.push([bundle, id]);

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    return true;
  }

  return getParents(global.parcelRequire, id).some(function (id) {
    return hmrAcceptCheck(global.parcelRequire, id);
  });
}

function hmrAcceptRun(bundle, id) {
  var cached = bundle.cache[id];
  bundle.hotData = {};

  if (cached) {
    cached.hot.data = bundle.hotData;
  }

  if (cached && cached.hot && cached.hot._disposeCallbacks.length) {
    cached.hot._disposeCallbacks.forEach(function (cb) {
      cb(bundle.hotData);
    });
  }

  delete bundle.cache[id];
  bundle(id);
  cached = bundle.cache[id];

  if (cached && cached.hot && cached.hot._acceptCallbacks.length) {
    cached.hot._acceptCallbacks.forEach(function (cb) {
      cb();
    });

    return true;
  }
}
},{}]},{},["node_modules/parcel/src/builtins/hmr-runtime.js","Main.ts"], null)
//# sourceMappingURL=/Main.d562fc5b.js.map