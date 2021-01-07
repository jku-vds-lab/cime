/******/ (function(modules) { // webpackBootstrap
/******/ 	// The module cache
/******/ 	var installedModules = {};
/******/
/******/ 	// The require function
/******/ 	function __webpack_require__(moduleId) {
/******/
/******/ 		// Check if module is in cache
/******/ 		if(installedModules[moduleId]) {
/******/ 			return installedModules[moduleId].exports;
/******/ 		}
/******/ 		// Create a new module (and put it into the cache)
/******/ 		var module = installedModules[moduleId] = {
/******/ 			i: moduleId,
/******/ 			l: false,
/******/ 			exports: {}
/******/ 		};
/******/
/******/ 		// Execute the module function
/******/ 		modules[moduleId].call(module.exports, module, module.exports, __webpack_require__);
/******/
/******/ 		// Flag the module as loaded
/******/ 		module.l = true;
/******/
/******/ 		// Return the exports of the module
/******/ 		return module.exports;
/******/ 	}
/******/
/******/
/******/ 	// expose the modules object (__webpack_modules__)
/******/ 	__webpack_require__.m = modules;
/******/
/******/ 	// expose the module cache
/******/ 	__webpack_require__.c = installedModules;
/******/
/******/ 	// define getter function for harmony exports
/******/ 	__webpack_require__.d = function(exports, name, getter) {
/******/ 		if(!__webpack_require__.o(exports, name)) {
/******/ 			Object.defineProperty(exports, name, { enumerable: true, get: getter });
/******/ 		}
/******/ 	};
/******/
/******/ 	// define __esModule on exports
/******/ 	__webpack_require__.r = function(exports) {
/******/ 		if(typeof Symbol !== 'undefined' && Symbol.toStringTag) {
/******/ 			Object.defineProperty(exports, Symbol.toStringTag, { value: 'Module' });
/******/ 		}
/******/ 		Object.defineProperty(exports, '__esModule', { value: true });
/******/ 	};
/******/
/******/ 	// create a fake namespace object
/******/ 	// mode & 1: value is a module id, require it
/******/ 	// mode & 2: merge all properties of value into the ns
/******/ 	// mode & 4: return value when already ns object
/******/ 	// mode & 8|1: behave like require
/******/ 	__webpack_require__.t = function(value, mode) {
/******/ 		if(mode & 1) value = __webpack_require__(value);
/******/ 		if(mode & 8) return value;
/******/ 		if((mode & 4) && typeof value === 'object' && value && value.__esModule) return value;
/******/ 		var ns = Object.create(null);
/******/ 		__webpack_require__.r(ns);
/******/ 		Object.defineProperty(ns, 'default', { enumerable: true, value: value });
/******/ 		if(mode & 2 && typeof value != 'string') for(var key in value) __webpack_require__.d(ns, key, function(key) { return value[key]; }.bind(null, key));
/******/ 		return ns;
/******/ 	};
/******/
/******/ 	// getDefaultExport function for compatibility with non-harmony modules
/******/ 	__webpack_require__.n = function(module) {
/******/ 		var getter = module && module.__esModule ?
/******/ 			function getDefault() { return module['default']; } :
/******/ 			function getModuleExports() { return module; };
/******/ 		__webpack_require__.d(getter, 'a', getter);
/******/ 		return getter;
/******/ 	};
/******/
/******/ 	// Object.prototype.hasOwnProperty.call
/******/ 	__webpack_require__.o = function(object, property) { return Object.prototype.hasOwnProperty.call(object, property); };
/******/
/******/ 	// __webpack_public_path__
/******/ 	__webpack_require__.p = "";
/******/
/******/
/******/ 	// Load entry module and return exports
/******/ 	return __webpack_require__(__webpack_require__.s = "./src/components/workers/worker_healthcheck.tsx");
/******/ })
/************************************************************************/
/******/ ({

/***/ "./src/components/workers/worker_healthcheck.tsx":
/*!*******************************************************!*\
  !*** ./src/components/workers/worker_healthcheck.tsx ***!
  \*******************************************************/
/*! no static exports found */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


Object.defineProperty(exports, "__esModule", {
  value: true
});

var backend_utils = __webpack_require__(/*! ../../utils/backend-connect */ "./src/utils/backend-connect.ts");

var backendRunning = false;

function check() {
  fetch(backend_utils.BASE_URL + '/healthcheck', {
    method: 'GET'
  }).then(function () {
    var context = self;

    if (backendRunning == false) {
      context.postMessage(true);
    }

    backendRunning = true;
  }).catch(function () {
    var context = self;

    if (backendRunning == true) {
      context.postMessage(false);
    }

    backendRunning = false;
  });
}

check();

/***/ }),

/***/ "./src/utils/backend-connect.ts":
/*!**************************************!*\
  !*** ./src/utils/backend-connect.ts ***!
  \**************************************/
/*! no static exports found */
/***/ (function(module, exports, __webpack_require__) {

"use strict";


var __awaiter = this && this.__awaiter || function (thisArg, _arguments, P, generator) {
  function adopt(value) {
    return value instanceof P ? value : new P(function (resolve) {
      resolve(value);
    });
  }

  return new (P || (P = Promise))(function (resolve, reject) {
    function fulfilled(value) {
      try {
        step(generator.next(value));
      } catch (e) {
        reject(e);
      }
    }

    function rejected(value) {
      try {
        step(generator["throw"](value));
      } catch (e) {
        reject(e);
      }
    }

    function step(result) {
      result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected);
    }

    step((generator = generator.apply(thisArg, _arguments || [])).next());
  });
};

Object.defineProperty(exports, "__esModule", {
  value: true
});
exports.calculate_hdbscan_clusters = exports.upload_sdf_file = exports.get_mcs_from_smiles_list = exports.get_structures_from_smiles_list = exports.get_structure_from_smiles = exports.BASE_URL = void 0;
exports.BASE_URL = 'http://127.0.0.1:8080'; // for local
// export const BASE_URL = ''; // for AWS
// export const BASE_URL = 'https://chemvis.caleydoapp.org'; // for netlify
// export const BASE_URL = 'http://127.0.0.1:5000';
// export const BASE_URL = 'http://caleydoapp.org:32819';

var smiles_cache = {};

function get_smiles_cache(smiles) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee() {
    return regeneratorRuntime.wrap(function _callee$(_context) {
      while (1) {
        switch (_context.prev = _context.next) {
          case 0:
            return _context.abrupt("return", smiles_cache[smiles]);

          case 1:
          case "end":
            return _context.stop();
        }
      }
    }, _callee);
  }));
}

function get_structure_from_smiles(smiles) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee2() {
    var formData;
    return regeneratorRuntime.wrap(function _callee2$(_context2) {
      while (1) {
        switch (_context2.prev = _context2.next) {
          case 0:
            if (!smiles_cache[smiles]) {
              _context2.next = 2;
              break;
            }

            return _context2.abrupt("return", get_smiles_cache(smiles));

          case 2:
            formData = new FormData();
            formData.append('smiles', smiles);
            return _context2.abrupt("return", fetch(exports.BASE_URL + '/get_mol_img', {
              method: 'POST',
              body: formData
            }).then(function (response) {
              return response.text();
            }).then(function (data) {
              smiles_cache[smiles] = data;
              return data;
            }).catch(function (error) {
              console.error(error);
            }));

          case 5:
          case "end":
            return _context2.stop();
        }
      }
    }, _callee2);
  }));
}

exports.get_structure_from_smiles = get_structure_from_smiles;

function get_structures_from_smiles_list(formData) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee3() {
    return regeneratorRuntime.wrap(function _callee3$(_context3) {
      while (1) {
        switch (_context3.prev = _context3.next) {
          case 0:
            return _context3.abrupt("return", fetch(exports.BASE_URL + '/get_mol_imgs', {
              method: 'POST',
              body: formData
            }).then(function (response) {
              return response.json();
            }).catch(function (error) {
              console.error(error);
            }));

          case 1:
          case "end":
            return _context3.stop();
        }
      }
    }, _callee3);
  }));
}

exports.get_structures_from_smiles_list = get_structures_from_smiles_list;

function get_mcs_from_smiles_list(formData) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee4() {
    return regeneratorRuntime.wrap(function _callee4$(_context4) {
      while (1) {
        switch (_context4.prev = _context4.next) {
          case 0:
            return _context4.abrupt("return", fetch(exports.BASE_URL + '/get_common_mol_img', {
              method: 'POST',
              body: formData
            }).then(function (response) {
              return response.text();
            }).catch(function (error) {
              console.error(error);
            }));

          case 1:
          case "end":
            return _context4.stop();
        }
      }
    }, _callee4);
  }));
}

exports.get_mcs_from_smiles_list = get_mcs_from_smiles_list;

function upload_sdf_file(file) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee5() {
    var formData_file;
    return regeneratorRuntime.wrap(function _callee5$(_context5) {
      while (1) {
        switch (_context5.prev = _context5.next) {
          case 0:
            // upload the sdf file to the server
            // the response is a unique filename that can be used to make further requests
            formData_file = new FormData();
            formData_file.append('myFile', file);
            return _context5.abrupt("return", fetch(exports.BASE_URL + '/upload_sdf', {
              method: 'POST',
              body: formData_file
            }).then(function (response) {
              return response.json();
            }).catch(function (error) {
              console.error(error);
            }));

          case 3:
          case "end":
            return _context5.stop();
        }
      }
    }, _callee5);
  }));
}

exports.upload_sdf_file = upload_sdf_file;

function calculate_hdbscan_clusters(X) {
  return __awaiter(this, void 0, void 0, /*#__PURE__*/regeneratorRuntime.mark(function _callee6() {
    return regeneratorRuntime.wrap(function _callee6$(_context6) {
      while (1) {
        switch (_context6.prev = _context6.next) {
          case 0:
            return _context6.abrupt("return", fetch(exports.BASE_URL + '/segmentation', {
              method: 'POST',
              body: JSON.stringify(X)
            }).then(function (response) {
              return response.json();
            }));

          case 1:
          case "end":
            return _context6.stop();
        }
      }
    }, _callee6);
  }));
}

exports.calculate_hdbscan_clusters = calculate_hdbscan_clusters;

/***/ })

/******/ });
//# sourceMappingURL=healthcheck.js.map