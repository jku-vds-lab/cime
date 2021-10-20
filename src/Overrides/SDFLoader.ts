import { CIMEBackendFromEnv } from "../Backend/CIMEBackend";
import { trackPromise } from "react-promise-tracker";
import {
  AVector,
  CSVLoader,
  Dataset,
  DatasetType,
  IVector,
  Loader,
  useCancellablePromise,
  // DatasetEntry,
} from "projection-space-explorer";

var d3v5 = require("d3");

function convertFromCSV(vectors) {
  return vectors.map((vector) => {
    return AVector.create(vector);
  });
}

export class SDFLoader implements Loader {
  vectors: IVector[] = [];
  datasetType: DatasetType = DatasetType.None;

  loading_area = "global_loading_indicator";

  resolvePath(
    entry: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    abort_controller?
  ) {
    if (entry.uploaded) {
      // TODO: Instead of localstorage store it in state?
      // localStorage.setItem("id", entry.path);
      // use file that is already uploaded to backend
      this.loadCSV(
        finished,
        entry,
        cancellablePromise,
        modifiers,
        abort_controller
      );
    } else {
      trackPromise(
        fetch(entry.path, { signal: abort_controller?.signal })
          .then((response) => response.blob())
          .then((result) =>
            this.resolveContent(
              result,
              finished,
              cancellablePromise,
              modifiers,
              abort_controller
            )
          )
          .catch((error) => {
            console.log(error);
          }),
        this.loading_area
      );
    }
  }

  resolveContent(
    file: any,
    finished: (dataset: Dataset) => void,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    controller?
  ) {
    const promise = cancellablePromise
      ? cancellablePromise(
          CIMEBackendFromEnv.upload_sdf_file(file, controller),
          controller
        )
      : CIMEBackendFromEnv.upload_sdf_file(file, controller);
    trackPromise(
      promise
        .then((uploaded) => {
          console.log("UPLOADED", uploaded);
          this.loadCSV(
            finished,
            { display: "", type: this.datasetType, path: uploaded.id },
            cancellablePromise,
            modifiers,
            controller
          );
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loading_area
    );
  }

  loadCSV(
    finished: (dataset: Dataset) => void,
    entry,
    cancellablePromise?: ReturnType<
      typeof useCancellablePromise
    >["cancellablePromise"],
    modifiers?: string,
    controller?
  ) {
    // request the server to return a csv file using the unique filename
    const path =
      CIMEBackendFromEnv.baseUrl + "/get_csv/" + entry.path + "/" + modifiers;
    const promise = cancellablePromise
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
        });
    trackPromise(
      promise
        .then((vectors) => {
          this.vectors = convertFromCSV(vectors);
          this.datasetType = DatasetType.Chem;
          new CSVLoader().resolve(
            finished,
            this.vectors,
            this.datasetType,
            entry
          );
        })
        .catch((error) => {
          console.log(error);
        }),
      this.loading_area
    );
  }
}
