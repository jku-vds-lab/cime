var __assign = (this && this.__assign) || function () {
    __assign = Object.assign || function(t) {
        for (var s, i = 1, n = arguments.length; i < n; i++) {
            s = arguments[i];
            for (var p in s) if (Object.prototype.hasOwnProperty.call(s, p))
                t[p] = s[p];
        }
        return t;
    };
    return __assign.apply(this, arguments);
};
var __awaiter = (this && this.__awaiter) || function (thisArg, _arguments, P, generator) {
    function adopt(value) { return value instanceof P ? value : new P(function (resolve) { resolve(value); }); }
    return new (P || (P = Promise))(function (resolve, reject) {
        function fulfilled(value) { try { step(generator.next(value)); } catch (e) { reject(e); } }
        function rejected(value) { try { step(generator["throw"](value)); } catch (e) { reject(e); } }
        function step(result) { result.done ? resolve(result.value) : adopt(result.value).then(fulfilled, rejected); }
        step((generator = generator.apply(thisArg, _arguments || [])).next());
    });
};
var __generator = (this && this.__generator) || function (thisArg, body) {
    var _ = { label: 0, sent: function() { if (t[0] & 1) throw t[1]; return t[1]; }, trys: [], ops: [] }, f, y, t, g;
    return g = { next: verb(0), "throw": verb(1), "return": verb(2) }, typeof Symbol === "function" && (g[Symbol.iterator] = function() { return this; }), g;
    function verb(n) { return function (v) { return step([n, v]); }; }
    function step(op) {
        if (f) throw new TypeError("Generator is already executing.");
        while (_) try {
            if (f = 1, y && (t = op[0] & 2 ? y["return"] : op[0] ? y["throw"] || ((t = y["return"]) && t.call(y), 0) : y.next) && !(t = t.call(y, op[1])).done) return t;
            if (y = 0, t) op = [op[0] & 2, t.value];
            switch (op[0]) {
                case 0: case 1: t = op; break;
                case 4: _.label++; return { value: op[1], done: false };
                case 5: _.label++; y = op[1]; op = [0]; continue;
                case 7: op = _.ops.pop(); _.trys.pop(); continue;
                default:
                    if (!(t = _.trys, t = t.length > 0 && t[t.length - 1]) && (op[0] === 6 || op[0] === 2)) { _ = 0; continue; }
                    if (op[0] === 3 && (!t || (op[1] > t[0] && op[1] < t[3]))) { _.label = op[1]; break; }
                    if (op[0] === 6 && _.label < t[1]) { _.label = t[1]; t = op; break; }
                    if (t && _.label < t[2]) { _.label = t[2]; _.ops.push(op); break; }
                    if (t[2]) _.ops.pop();
                    _.trys.pop(); continue;
            }
            op = body.call(thisArg, _);
        } catch (e) { op = [6, e]; y = 0; } finally { f = t = 0; }
        if (op[0] & 5) throw op[1]; return { value: op[0] ? op[1] : void 0, done: true };
    }
};
var CIMEBackend = /** @class */ (function () {
    function CIMEBackend(baseUrl, fetchParams) {
        var _this = this;
        if (fetchParams === void 0) { fetchParams = {}; }
        this.baseUrl = baseUrl;
        this.fetchParams = fetchParams;
        this.smiles_cache = {};
        this.smiles_highlight_cache = {};
        this.cache = {};
        this.handleSmilesCache = function (smiles, highlight) {
            if (highlight === void 0) { highlight = false; }
            //already downloaded this image -> saved in smiles cache
            if (highlight) {
                return _this.smiles_highlight_cache[smiles];
            }
            else {
                return _this.smiles_cache[smiles];
            }
        };
        this.setSmilesCache = function (smiles, highlight, data) {
            if (highlight === void 0) { highlight = false; }
            if (highlight)
                _this.smiles_highlight_cache[smiles] = data;
            else
                _this.smiles_cache[smiles] = data;
        };
        this.async_cache = function (cached_data) { return __awaiter(_this, void 0, void 0, function () {
            return __generator(this, function (_a) {
                return [2 /*return*/, cached_data];
            });
        }); };
        this.handleCache = function (key) {
            if (_this.cache[key])
                return Object.assign(_this.cache[key]); // return copy of cached object
            return null;
        };
        this.setCache = function (key, value) {
            _this.cache[key] = value;
        };
        this.handleErrors = function (response) {
            if (!response.ok) {
                throw Error(response.statusText);
            }
            return response;
        };
        this.handleJSONErrors = function (data) {
            if (Object.keys(data).includes("error")) {
                alert(data["error"]);
            }
            return data;
        };
        this.deleteFile = function (filename) { return __awaiter(_this, void 0, void 0, function () {
            var path;
            return __generator(this, function (_a) {
                path = this.baseUrl + "/delete_file/" + filename;
                return [2 /*return*/, fetch(path, __assign(__assign({}, this.fetchParams), { method: "GET" }))
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .catch(function (error) {
                        alert("file could not be deleted. please, try again");
                        console.log(error);
                    })];
            });
        }); };
        this.getFiles = function () { return __awaiter(_this, void 0, void 0, function () {
            var path;
            return __generator(this, function (_a) {
                path = this.baseUrl + "/get_uploaded_files_list";
                return [2 /*return*/, fetch(path, __assign(__assign({}, this.fetchParams), { method: "GET" }))
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .catch(function (error) {
                        // alert("could not load uploaded filenames.")
                        console.log(error);
                    })];
            });
        }); };
        this.getDifferenceHighlight = function (smilesA, smilesB, controller) { return __awaiter(_this, void 0, void 0, function () {
            var formData, path, my_fetch;
            return __generator(this, function (_a) {
                formData = new FormData();
                formData.append("smilesA", smilesA);
                formData.append("smilesB", smilesB);
                path = this.baseUrl + "/get_difference_highlight";
                if (controller) {
                    my_fetch = fetch(path, __assign(__assign({}, this.fetchParams), { method: "POST", body: formData, signal: controller.signal }));
                }
                else {
                    my_fetch = fetch(path, __assign(__assign({}, this.fetchParams), { method: "POST", body: formData }));
                }
                return [2 /*return*/, my_fetch
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .then(function (data) {
                        console.log(data);
                        return data["data"];
                    })
                        .catch(function (error) {
                        // alert("could not load structure");
                        console.log(error);
                    })];
            });
        }); };
        this.getStructureFromSmiles = function (id, smiles, highlight, controller) {
            if (highlight === void 0) { highlight = false; }
            var cached_data = _this.handleSmilesCache(smiles, highlight);
            if (cached_data) {
                return _this.async_cache(cached_data);
            }
            var formData = new FormData();
            formData.append("smiles", smiles);
            formData.append("filename", id === null || id === void 0 ? void 0 : id.toString());
            var path = _this.baseUrl + "/get_mol_img";
            if (highlight) {
                path += "/highlight";
            }
            return fetch(path, __assign(__assign({}, _this.fetchParams), { method: "POST", body: formData, signal: controller === null || controller === void 0 ? void 0 : controller.signal }))
                .then(_this.handleErrors)
                .then(function (response) { return response.json(); })
                .then(_this.handleJSONErrors)
                .then(function (data) {
                _this.setSmilesCache(smiles, highlight, data["data"]);
                return data["data"];
            })
                .catch(function (error) {
                // alert("could not load structure");
                console.log(error);
            });
        };
        this.getStructuresFromSmilesList = function (formData, controller) { return __awaiter(_this, void 0, void 0, function () {
            return __generator(this, function (_a) {
                return [2 /*return*/, fetch(this.baseUrl + "/get_mol_imgs", __assign(__assign({}, this.fetchParams), { method: "POST", body: formData, signal: controller === null || controller === void 0 ? void 0 : controller.signal }))
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .then(function (data) {
                        if (data["error_smiles"].length > 0) {
                            alert("following smiles couldn not be parsed: " + data["error_smiles"]);
                        }
                        return data;
                    })
                        .catch(function (error) {
                        if (error.name === "AbortError") {
                            console.log("Fetch aborted");
                        }
                        else {
                            alert("could not load structures");
                            console.log(error);
                        }
                    })];
            });
        }); };
        this.getMCSFromSmilesList = function (formData, controller) { return __awaiter(_this, void 0, void 0, function () {
            var my_fetch;
            return __generator(this, function (_a) {
                if (controller) {
                    my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
                        method: "POST",
                        body: formData,
                        signal: controller === null || controller === void 0 ? void 0 : controller.signal,
                    });
                }
                else {
                    my_fetch = fetch(this.baseUrl + "/get_common_mol_img", {
                        method: "POST",
                        body: formData,
                    });
                }
                return [2 /*return*/, my_fetch
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .then(function (response) { return response["data"]; })
                        .catch(function (error) {
                        // alert("could not get maximum common substructure")
                        console.log(error);
                    })];
            });
        }); };
        this.getSubstructureCount = function (smiles_list, filter) { return __awaiter(_this, void 0, void 0, function () {
            var formData;
            return __generator(this, function (_a) {
                formData = new FormData();
                formData.append("smiles_list", smiles_list);
                formData.append("filter_smiles", filter);
                return [2 /*return*/, fetch(this.baseUrl + "/get_substructure_count", {
                        method: "POST",
                        body: formData,
                    })
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .then(function (data) {
                        if (Object.keys(data).includes("substructure_counts"))
                            return data["substructure_counts"];
                        else
                            throw Error("Backend responded with error: " + data["error"]);
                    })
                        .catch(function (error) {
                        alert("could not find substructure match");
                        console.log(error);
                    })];
            });
        }); };
        this.upload_sdf_file = function (file, controller) { return __awaiter(_this, void 0, void 0, function () {
            var formData_file, promise;
            return __generator(this, function (_a) {
                formData_file = new FormData();
                formData_file.append("myFile", file);
                formData_file.append("file_size", file.size);
                promise = fetch(this.baseUrl + "/upload_sdf", __assign(__assign({}, this.fetchParams), { method: "POST", body: formData_file, signal: controller === null || controller === void 0 ? void 0 : controller.signal }))
                    .then(this.handleErrors)
                    .then(function (response) { return response.json(); })
                    .then(this.handleJSONErrors)
                    .catch(function (error) {
                    if (error.name === "AbortError") {
                        console.log("Fetch aborted");
                    }
                    else {
                        alert("error when uploading file. it might be too big");
                        console.log(error);
                    }
                });
                return [2 /*return*/, promise];
            });
        }); };
        this.getRepresentationList = function (refresh, id, controller) {
            if (refresh === void 0) { refresh = false; }
            return __awaiter(_this, void 0, void 0, function () {
                var cached_data, path;
                var _this = this;
                return __generator(this, function (_a) {
                    if (!refresh) {
                        cached_data = this.handleCache("representation_list_" + id);
                        if (cached_data && cached_data["rep_list"].length > 0) {
                            return [2 /*return*/, this.async_cache(cached_data)];
                        }
                    }
                    path = this.baseUrl + ("/get_atom_rep_list/" + id);
                    return [2 /*return*/, fetch(path, __assign(__assign({}, this.fetchParams), { method: "GET", signal: controller === null || controller === void 0 ? void 0 : controller.signal }))
                            .then(this.handleErrors)
                            .then(function (response) { return response.json(); })
                            .then(this.handleJSONErrors)
                            .then(function (data) {
                            _this.setCache("representation_list_" + id, data);
                            return data;
                        })
                            .catch(function (error) {
                            // alert("error when loading representation list")
                            console.log(error);
                        })];
                });
            });
        };
        this.calculateHDBScanClusters = function (X, min_cluster_size, min_cluster_samples, allow_single_cluster) { return __awaiter(_this, void 0, void 0, function () {
            var formData;
            return __generator(this, function (_a) {
                formData = new FormData();
                formData.append("min_cluster_size", min_cluster_size);
                formData.append("min_cluster_samples", min_cluster_samples);
                formData.append("allow_single_cluster", allow_single_cluster);
                formData.append("X", X);
                return [2 /*return*/, fetch(this.baseUrl + "/segmentation", {
                        method: "POST",
                        body: formData,
                    })
                        .then(this.handleErrors)
                        .then(function (response) { return response.json(); })
                        .then(this.handleJSONErrors)
                        .catch(function (error) {
                        alert("error when calculating clusters");
                        console.log(error);
                    })];
            });
        }); };
    }
    return CIMEBackend;
}());
export { CIMEBackend };
// Use the environment variables defined in the .env file
if (!process.env.REACT_APP_CIME_BACKEND_URL) {
    console.error("The ENV-variable REACT_APP_CIME_BACKEND_URL must be set.");
}
export var CIMEBackendFromEnv = new CIMEBackend(process.env.REACT_APP_CIME_BACKEND_URL, {
// credentials: "omit",
});
//# sourceMappingURL=CIMEBackend.js.map