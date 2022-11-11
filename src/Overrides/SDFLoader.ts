import { CIMEBackendFromEnv } from '../Backend/CIMEBackend';
import { trackPromise } from 'react-promise-tracker';
import { AVector, CSVLoader, Dataset, DatasetType, DefaultFeatureLabel, IVector, Loader, useCancellablePromise, Profiler, VectView } from 'projection-space-explorer';
import * as Papa from 'papaparse';

export class SDFLoader implements Loader {
  vectors: IVector[] = [];
  datasetType: DatasetType = DatasetType.None;

  loadingArea = 'global_loading_indicator';

  resolvePath(
    entry: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    modifiers?: string,
    controller?: AbortController,
  ) {
    if (entry.uploaded) {
      // TODO: Instead of localstorage store it in state?
      // localStorage.setItem("id", entry.path);
      // use file that is already uploaded to backend
      this.loadCSV(finished, entry, cancellablePromise, modifiers, controller);
    } else {
      trackPromise(
        fetch(entry.path, { signal: controller?.signal })
          .then((response) => response.blob())
          .then((result) => this.resolveContent(result, finished, cancellablePromise, modifiers, controller))
          .catch((error) => {
            console.log(error);
          }),
        this.loadingArea,
      );
    }
  }

  resolveContent(
    file: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    modifiers?: string,
    controller?: AbortController,
  ) {
    const promise = cancellablePromise
      ? cancellablePromise(CIMEBackendFromEnv.upload_sdf_file(file, controller), controller)
      : CIMEBackendFromEnv.upload_sdf_file(file, controller);
    trackPromise(
      promise
        .then((uploaded) => {
          //console.log("UPLOADED", uploaded);
          this.loadCSV(finished, { display: '', type: this.datasetType, path: uploaded.id }, cancellablePromise, modifiers, controller);
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loadingArea,
    );
  }

  loadCSV(
    finished: (dataset: Dataset) => void,
    entry,
    cancellablePromise?: ReturnType<typeof useCancellablePromise>['cancellablePromise'],
    modifiers?: string,
    controller?: AbortController,
    onError?: (error: any) => void,
  ) {
    const profiler = new Profiler();
    // request the server to return a csv file using the unique filename
    const path = CIMEBackendFromEnv.baseUrl + '/get_csv/' + entry.path + '/' + modifiers;
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

    const metaInformation = {};
    const resultingVectors = [];

    const promise = new Promise((resolve) => {
      Papa.parse<IVector>(path, {
        download: true,
        header: true,
        skipEmptyLines: true,
        transformHeader: (value, o) => {
          const json = value.match(/[{].*[}]/);

          if (json != null) {
            const cutHeader = value.substring(0, value.length - json[0].length);
            metaInformation[cutHeader] = JSON.parse(json[0]);
            return cutHeader;
          } else {
            metaInformation[value] = { featureLabel: DefaultFeatureLabel };
            return value;
          }
        },
        downloadRequestHeaders: {
          credentials: 'same-origin',
        },
        step: function(results) {
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
    

    trackPromise(
      promise
        .then(() => {
          profiler.profile("loaded data");
          this.vectors = resultingVectors;

          profiler.profile('converted to vectors');
         
          this.datasetType = DatasetType.Chem;

          new CSVLoader().resolve(finished, this.vectors, this.datasetType, entry, metaInformation).catch((error) => {
            onError(error);
          });
        })
        .catch((error) => {
          if (onError) {
            onError(error);
          }
        }),
      this.loadingArea,
    );
  }
}
