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
import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import { trackPromise } from "react-promise-tracker";
import { AVector, CSVLoader, DatasetType,
// DatasetEntry,
 } from "projection-space-explorer";
import * as d3v5 from "d3v5";
function convertFromCSV(vectors) {
    return vectors.map(function (vector) {
        return AVector.create(vector);
    });
}
var SDFLoader = /** @class */ (function () {
    function SDFLoader() {
        this.vectors = [];
        this.datasetType = DatasetType.None;
        this.loadingArea = "global_loading_indicator";
    }
    SDFLoader.prototype.resolvePath = function (entry, finished, cancellablePromise, modifiers, controller) {
        var _this = this;
        if (entry.uploaded) {
            // TODO: Instead of localstorage store it in state?
            // localStorage.setItem("id", entry.path);
            // use file that is already uploaded to backend
            this.loadCSV(finished, entry, cancellablePromise, modifiers, controller);
        }
        else {
            trackPromise(fetch(entry.path, { signal: controller === null || controller === void 0 ? void 0 : controller.signal })
                .then(function (response) { return response.blob(); })
                .then(function (result) {
                return _this.resolveContent(result, finished, cancellablePromise, modifiers, controller);
            })
                .catch(function (error) {
                console.log(error);
            }), this.loadingArea);
        }
    };
    SDFLoader.prototype.resolveContent = function (file, finished, cancellablePromise, modifiers, controller) {
        var _this = this;
        var promise = cancellablePromise
            ? cancellablePromise(CIMEBackendFromEnv.upload_sdf_file(file, controller), controller)
            : CIMEBackendFromEnv.upload_sdf_file(file, controller);
        trackPromise(promise
            .then(function (uploaded) {
            //console.log("UPLOADED", uploaded);
            _this.loadCSV(finished, { display: "", type: _this.datasetType, path: uploaded.id }, cancellablePromise, modifiers, controller);
        })
            .catch(function (error) {
            console.log(error);
        }), this.loadingArea);
    };
    SDFLoader.prototype.loadCSV = function (finished, entry, cancellablePromise, modifiers, controller) {
        var _this = this;
        // request the server to return a csv file using the unique filename
        var path = CIMEBackendFromEnv.baseUrl + "/get_csv/" + entry.path + "/" + modifiers;
        var promise = cancellablePromise
            ? cancellablePromise(d3v5.csv(path, __assign(__assign({}, CIMEBackendFromEnv.fetchParams), { signal: controller === null || controller === void 0 ? void 0 : controller.signal })), controller)
            : d3v5.csv(path, __assign(__assign({}, CIMEBackendFromEnv.fetchParams), { signal: controller === null || controller === void 0 ? void 0 : controller.signal }));
        trackPromise(promise
            .then(function (vectors) {
            _this.vectors = convertFromCSV(vectors);
            _this.datasetType = DatasetType.Chem;
            new CSVLoader().resolve(finished, _this.vectors, _this.datasetType, entry);
        })
            .catch(function (error) {
            console.log(error);
        }), this.loadingArea);
    };
    return SDFLoader;
}());
export { SDFLoader };
//# sourceMappingURL=SDFLoader.js.map