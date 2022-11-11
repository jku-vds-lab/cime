import { CIMEBackendFromEnv } from '../Backend/CIMEBackend';
import { trackPromise } from 'react-promise-tracker';
import { CSVLoader, DatasetType, DefaultFeatureLabel, Profiler, VectView } from 'projection-space-explorer';
import * as Papa from 'papaparse';
var SDFLoader = /** @class */ (function () {
    function SDFLoader() {
        this.vectors = [];
        this.datasetType = DatasetType.None;
        this.loadingArea = 'global_loading_indicator';
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
                .then(function (result) { return _this.resolveContent(result, finished, cancellablePromise, modifiers, controller); })
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
            _this.loadCSV(finished, { display: '', type: _this.datasetType, path: uploaded.id }, cancellablePromise, modifiers, controller);
        })
            .catch(function (error) {
            console.log(error);
        }), this.loadingArea);
    };
    SDFLoader.prototype.loadCSV = function (finished, entry, cancellablePromise, modifiers, controller, onError) {
        var _this = this;
        var profiler = new Profiler();
        // request the server to return a csv file using the unique filename
        var path = CIMEBackendFromEnv.baseUrl + '/get_csv/' + entry.path + '/' + modifiers;
        /**const promise = cancellablePromise
          ? cancellablePromise(
              d3v5.csv(path, {
                ...CIMEBackendFromEnv.fetchParams,
                signal: controller?.signal,
              }),
              controller
            )
          : d3v5.csv(path, {
              ...CIMEBackendFromEnv.fetchParams,
              signal: controller?.signal,
            });**/
        profiler.profile('start');
        var metaInformation = {};
        var resultingVectors = [];
        var promise = new Promise(function (resolve) {
            Papa.parse(path, {
                download: true,
                header: true,
                skipEmptyLines: true,
                transformHeader: function (value, o) {
                    var json = value.match(/[{].*[}]/);
                    if (json != null) {
                        var cutHeader = value.substring(0, value.length - json[0].length);
                        metaInformation[cutHeader] = JSON.parse(json[0]);
                        return cutHeader;
                    }
                    else {
                        metaInformation[value] = { featureLabel: DefaultFeatureLabel };
                        return value;
                    }
                },
                downloadRequestHeaders: {
                    credentials: 'same-origin',
                },
                step: function (results) {
                    // Inject default attributes of vectors
                    results.data.objectType = 'vector';
                    results.data.__meta__ = new VectView();
                    resultingVectors.push(results.data);
                },
                complete: function (results) {
                    resolve(resultingVectors);
                },
            });
        });
        trackPromise(promise
            .then(function () {
            profiler.profile("loaded data");
            _this.vectors = resultingVectors;
            profiler.profile('converted to vectors');
            _this.datasetType = DatasetType.Chem;
            new CSVLoader().resolve(finished, _this.vectors, _this.datasetType, entry, metaInformation).catch(function (error) {
                onError(error);
            });
        })
            .catch(function (error) {
            if (onError) {
                onError(error);
            }
        }), this.loadingArea);
    };
    return SDFLoader;
}());
export { SDFLoader };
//# sourceMappingURL=SDFLoader.js.map